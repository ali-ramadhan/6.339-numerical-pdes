function rho_xt = q1()
% Problem parameters.
rho_max = 1.0;
v_max = 1.0;
x_max = 10;

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

% Numerical scheme parameters
N = 500;
x = linspace(0,x_max,N+1);
dx = x(2:end) - x(1:end-1);
dt = 0.01;
n_step = 2/dt;

% initial condition
rho_0 = 0.2*rho_max*ones(N+1,1); % light traffic
% rho_0 = 0.8*rho_max*ones(N+1,1); % heavy traffic
rho = rho_0;

% We will store rho(x,t) for all cell centers x_i and 
rho_xt = zeros(N+1, n_step+1);
rho_xt(:,1) = rho_0;

for i=1:n_step
    % rho = zeros(N+1,1);
    rho(1) = rho_0(1); % Impose left BC: rho(0,t) = rho_0.
    
    F_imh = zeros(N+1,1); % F_{i-1/2}
    F_iph = zeros(N+1,1); % F_{i+1/2}
    for j=2:N
        F_imh(j) = F(rho(j-1), rho(j));
        F_iph(j) = F(rho(j), rho(j+1));
    end
    
    % Simulate accident at x=5 for t<1 by setting the flux going in and out
    % of cell N/2 to 0.
    if i*dt < 1
        F_imh(int64(N/2)+1) = 0;
        F_iph(int64(N/2)) = 0;
    end
    
    size(rho)
    size(dt./dx)
    size(F_iph - F_imh)
    % First-order finite volume scheme
    rho = rho - (dt./dx) .* (F_iph - F_imh);
    
    rho_xt(:,i+1) = rho;
end

end