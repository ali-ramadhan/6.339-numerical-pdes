% * k = kappa and possible options for flux_scheme are 'none', 'minmod',
% 'superbee', or 'vanLeer'.
% * rho_init must be between 0 and rho_max.
% * n is the number of lanes.
function [x, t, rho_xt] = q3a(rho_init, n, k, flux_scheme)
% Problem parameters.
rho_max = 1.0;
v_max = 1.0;
x_max = 10;
alpha = 0.1;
% n = 3; % number of lanes

% Numerical scheme parameters
N = 500; % Keep even so we can have a middle cell N/2 for the accident.
x = linspace(0,x_max,N+1);
dx = [x(2:end) - x(1:end-1), x_max/N];
dt = 0.01;
n_step = 2/dt;
t = linspace(0,n_step*dt,n_step+1);

% Flux reconstruction function
function F = F(rho_L, rho_R)
    % Flux function
    function f = f(rho)
        f = v_max * (rho - (rho^3 / rho_max^2));
    end

    % Reconstruct flux using Godunov's scheme with the exact solution to
    % the Riemann problem.
    if rho_L <= rho_R
        % Will treat the rho_L == rho_R equal case here. It occurs with
        % flat initial conditions at the first time step.
        if rho_R <= rho_max/sqrt(3)
            F = f(rho_L);
        elseif rho_L >= rho_max/sqrt(3)
            F = f(rho_R);
        elseif rho_L < rho_max/sqrt(3) && rho_max/sqrt(3) < rho_R
            F = min(f(rho_L), f(rho_R));
        end
    elseif rho_L > rho_R
        if rho_L <= rho_max/sqrt(3)
            F = f(rho_L);
        elseif rho_R >= rho_max/sqrt(3)
            F = f(rho_R);
        elseif rho_L > rho_max/sqrt(3) && rho_max/sqrt(3) > rho_R
            F = f(rho_max/sqrt(3));
        end
    end
end

% Initial condition setup
rho_0 = rho_init*ones(n*(N+1),1);

% Derivative of rho for ode45.
function drhodt = q2odefun(t, rho, rho_0, dx)
    % Flux limiter
    function Psi = Psi(R_i)
        if strcmp(flux_scheme, 'minmod')
            Psi = max(0, min(R_i, 1));
        elseif strcmp(flux_scheme, 'superbee')
            Psi = max([0, min(2*R_i, 1), min(R_i, 2)]);
        elseif strcmp(flux_scheme, 'vanLeer')
            Psi = (R_i + abs(R_i))/(1 + abs(R_i));
        end
    end
    
    % Convert rho from n*(N+1)-vector to (N+1)xn matrix.
    rho = reshape(rho, [N+1,n]);
    
    rho(1,:) = rho_0(1,:); % Impose left BC: rho(0,t) = rho_0.
        
    % Add a ghost point on either side of rho to allow us to reconstruct
    % the flux at every grid point.
    rho = [rho(1,:); rho; rho(end,:)];
    
    % Limited kappa-reconstruction interface values.
    rho_minus = zeros(N+3,n);
    rho_plus  = zeros(N+3,n);
    for l=1:n
        for j=2:N+1
            if abs(rho(j,l) - rho(j-1,l)) < 1e-6 || abs(rho(j+1,l) - rho(j,l)) < 1e-6
                % To avoid Inf when calculating 1/R_i. Plugging in a high value
                % like 1e+10 into any flux scheme will reproduce its asymptotic
                % behavior very closely.
                R_i = 1e-10;
            else
                R_i = (rho(j+1,l) - rho(j,l))/(rho(j,l) - rho(j-1,l));
            end

            rho_minus(j,l) = rho(j,l) + ((1-k)/4)*Psi(R_i)*(rho(j,l) - rho(j-1,l)) ...
                + ((1+k)/4)*Psi(1/R_i)*(rho(j+1,l) - rho(j,l));
            rho_plus(j,l) = rho(j,l) - ((1-k)/4)*Psi(1/R_i)*(rho(j+1,l) - rho(j,l)) ...
                - ((1+k)/4)*Psi(R_i)*(rho(j,l) - rho(j-1,l));
        end
    end

    F_imh = zeros(N+3,n); % F_{i-1/2}
    F_iph = zeros(N+3,n); % F_{i+1/2}
    for l=1:n
        for j=2:N
            F_imh(j,l) = F(rho_plus(j-1,l), rho_minus(j,l));
            F_iph(j,l) = F(rho_plus(j,l), rho_minus(j+1,l));
        end
    end
    
    % Discard ghost point values so we end up with a (N+1)-vector again.
    F_iph = F_iph(2:end-1,:);
    F_imh = F_imh(2:end-1,:);
    rho = rho(2:end-1,:);
    
    % Construct s matrix
    s = zeros(N+1,n);
    if n==1
        ; % No lane switching if there's only one lane.
    else
        for l=1:n
            if l==1
                s(:,l) = alpha*(rho(:,2) - rho(:,1));
            elseif l==n
                s(:,l) = alpha*(rho(:,n-1) - rho(:,n));
            else
                s(:,l) = alpha*(rho(:,l+1) - 2*rho(:,l) + rho(:,l-1));
            end
        end
    end
    
    % Construct the derivative vector d(rho)/dt.
    drhodt = zeros(N+1,n);
    for l=1:n
        drhodt(:,l) = (1./dx') .* (F_imh(:,l) - F_iph(:,l)) + s(:,l);
    end
    drhodt(1,:) = 0; % Impose left BC: rho(0,t) = rho_0.
    
    % Represent braking in the middle lane by an impulsive source term for
    % t < 0.01.
    if t < 0.01
        if n == 1 || n == 2
            l = 1; % Choose first lane if we only have 1 or 2.
        else
            l = int64(n/2); % Choose the middle lane (rounded down).
        end
        
        % You really want a Dirac delta function here but a large number
        % should do the trick. Well, apparently not too large otherwise rho
        % blows up above rho_max...
        S = 50;
        drhodt(int64(N/2),l) = drhodt(int64(N/2),l) + S;
        drhodt(int64(N/2)+1,l) = drhodt(int64(N/2)+1,l) - S;
    end
    
    % Convert d(rho)/dt matrix into a column vector.
    drhodt = reshape(drhodt, [n*(N+1),1]);
end

[t,y] = ode45(@(t,rho) q2odefun(t, rho, rho_0, dx), [0, n_step*dt], rho_0');
rho_xt = y';

end