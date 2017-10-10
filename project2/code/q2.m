function [x, t, rho_xt] = q1()
% Problem parameters.
rho_max = 1.0;
v_max = 1.0;
x_max = 10;

% Numerical scheme parameters
N = 500; % Keep even so we can have a middle cell N/2 for the accident.
x = linspace(0,x_max,N+1);
dx = [x(2:end) - x(1:end-1), x_max/N];
dt = 0.0001;
n_step = 2/dt;
t = linspace(0,n_step*dt,n_step+1);
k = -1; % kappa

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

% initial condition
rho_0 = 0.2*rho_max*ones(N+1,1); % light traffic
% rho_0 = 0.8*rho_max*ones(N+1,1); % heavy traffic
rho = rho_0;

% We will store rho(x,t) for all cell centers x_i and 
rho_xt = zeros(N+1, n_step+1);
rho_xt(:,1) = rho_0;

for i=1:n_step
    rho(1) = rho_0(1); % Impose left BC: rho(0,t) = rho_0.
    
    plot(x,rho);
    title(i*dt);
    drawnow;
    
    % Add ghost points
    rho = [rho(1); rho; rho(end)];
    
    F_imh = zeros(N+1,1); % F_{i-1/2}
    F_iph = zeros(N+1,1); % F_{i+1/2}
%     k1 = zeros(N+1,1);
%     k2 = zeros(N+1,1);
%     k3 = zeros(N+1,1);
%     k4 = zeros(N+1,1);
    for j=3:N+1
        % Limited \kappa-reconstruction interface values
        rho_imh_minus = rho(j) + ((1-k)/4)*(rho(j-1) - rho(j-2)) ...
            + ((1+k)/4)*(rho(j) - rho(j-1));
        rho_imh_plus = rho(j) - ((1-k)/4)*(rho(j+1) - rho(j)) ...
            - ((1+k)/4)*(rho(j) - rho(j-1));
        
        rho_iph_minus = rho(j) + ((1-k)/4)*(rho(j) - rho(j-1)) ...
            + ((1+k)/4)*(rho(j+1) - rho(j));
        rho_iph_plus = rho(j) - ((1-k)/4)*(rho(j+2) - rho(j+1)) ...
            - ((1+k)/4)*(rho(j+1) - rho(j));
        
        F_imh(j) = F(rho_imh_minus, rho_imh_plus);
        F_iph(j) = F(rho_iph_minus, rho_iph_plus);
        
        % Simulate accident at x=5 for t<1 by setting the flux going in and out
        % of cell N/2 to 0.
        if i*dt < 1
            F_imh(int64(N/2)+1) = 0;
            F_iph(int64(N/2)) = 0;
        end
        
%         k1(j) = (1/dx(j)) * (F_imh(j) - F_iph(j));
%         k2(j) = (1/dx(j)) * (f(rho_imh_plus + (dx(j)/2)*k1(j)) ...
%             - f(rho_iph_minus + (dx(j)/2)*k1(j)));
%         k3(j) = (1/dx(j)) * (f(rho_imh_plus + (dx(j)/2)*k2(j)) ...
%             - f(rho_iph_minus + (dx(j)/2)*k2(j)));
%         k4(j) = (1/dx(j)) * (f(rho_imh_plus + dx(j)*k3(j)) ...
%             - f(rho_iph_minus + dx(j)*k3(j)));
    end
    
    % Fourth-order Runge-Kutta method
    % rho = rho - (dt/6).*(k1 + 2*k2 + 2*k3 + k4);
    
    % Remove ghost points
    rho = rho(2:end-1);
    
    rho = ode45();
    
    rho = rho - (dt./dx)' .* (F_iph - F_imh);  
    rho_xt(:,i+1) = rho;
end

end