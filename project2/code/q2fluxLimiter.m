function [x, t, rho_xt] = q2fluxLimiter(k, flux_scheme)
% Problem parameters.
rho_max = 1.0;
v_max = 1.0;
x_max = 10;

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
% rho_0 = 0.2*rho_max*ones(N+1,1); % light traffic
rho_0 = 0.8*rho_max*ones(N+1,1); % heavy traffic
rho = rho_0;

% We will store rho(x,t) for all cell centers x_i and time steps t_n.
rho_xt = zeros(N+1, n_step+1);
rho_xt(:,1) = rho_0;

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
    
    rho(1) = rho_0(1); % Impose left BC: rho(0,t) = rho_0.
        
    % Add a ghost point on either side of rho to allow us to reconstruct
    % the flux at every grid point.
    rho = [rho(1); rho; rho(end)];
    
    % Limited kappa-reconstruction interface values.
    rho_minus = zeros(N+3,1);
    rho_plus  = zeros(N+3,1);
    for j=2:N+1
        R_i = (rho(j+1) - rho(j))/(rho(j) - rho(j-1));
        rho_minus(j) = rho(j) + ((1-k)/4)*Psi(R_i)*(rho(j) - rho(j-1)) ...
            + ((1+k)/4)*Psi(1/R_i)*(rho(j+1) - rho(j));
        rho_plus(j) = rho(j) - ((1-k)/4)*Psi(1/R_i)*(rho(j+1) - rho(j)) ...
            - ((1+k)/4)*Psi(R_i)*(rho(j) - rho(j-1));
    end

    F_imh = zeros(N+3,1); % F_{i-1/2}
    F_iph = zeros(N+3,1); % F_{i+1/2}
    for j=2:N
        F_imh(j) = F(rho_plus(j-1), rho_minus(j));
        F_iph(j) = F(rho_plus(j), rho_minus(j+1));
    end
    
%     rho = rho(2:end);
%     rho_minus = rho_minus(2:end);
%     rho_plus = rho_plus(2:end);
    % Discard ghost point values so we end up with a (N+1)-vector again.
    F_iph = F_iph(2:end-1);
    F_imh = F_imh(2:end-1);        
  
    % Simulate accident at x=5 for t<1 by setting the flux going in and out
    % of cell N/2 to 0.
    if t < 1
        F_imh(int64(N/2)+1) = 0;
        F_iph(int64(N/2)) = 0;
    end
    
    % Construct the derivative vector d(rho)/dt.
    drhodt = (1./dx') .* (F_imh - F_iph);
    drhodt(1) = 0; % Impose left BC: rho(0,t) = rho_0.
end

[t,y] = ode45(@(t,rho) q2odefun(t, rho, rho_0, dx), [0, n_step*dt], rho_0');
rho_xt = y';

%     k1 = zeros(N+1,1);
%     k2 = zeros(N+1,1);
%     k3 = zeros(N+1,1);
%     k4 = zeros(N+1,1);

%         k1(j) = (1/dx(j)) * (F_imh(j) - F_iph(j));
%         k2(j) = (1/dx(j)) * (f(rho_imh_plus + (dx(j)/2)*k1(j)) ...
%             - f(rho_iph_minus + (dx(j)/2)*k1(j)));
%         k3(j) = (1/dx(j)) * (f(rho_imh_plus + (dx(j)/2)*k2(j)) ...
%             - f(rho_iph_minus + (dx(j)/2)*k2(j)));
%         k4(j) = (1/dx(j)) * (f(rho_imh_plus + dx(j)*k3(j)) ...
%             - f(rho_iph_minus + dx(j)*k3(j)));

% Fourth-order Runge-Kutta method
% rho = rho - (dt/6).*(k1 + 2*k2 + 2*k3 + k4);

end