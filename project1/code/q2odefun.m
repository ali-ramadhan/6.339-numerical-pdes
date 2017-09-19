function dydt = q2odefun(t, y, M, F)
global N dx dy

rho_0 = 1;
c_0 = 1;
u_0 = sqrt(2);
L = 6;
H = 2;

% Reshape y vector into p',q' matrices.
pp = reshape(y(1:(N+1)*(M+1)), N+1, M+1);
qp = reshape(y((N+1)*(M+1)+1:end), [N+1, M+1]);

% Imposing boundary conditions at each wall.
for j=1:M+1
    pp(1,j) = 0; % p=p_0 (p'=0) at x=0
    pp(N+1,j) = pp(N,j); % dp'/dx = 0 at x=L

    qp(1,j) = 0; % q'=0 at x=0
    qp(N+1,j) = qp(N,j); % dq'/dx 0 at x=L
end

for i=1:N+1
    pp(i,M+1) = pp(i,M); % dp'/dy = 0 at y=H

    if i*dx < L/4 || i*dx > L/2
        pp(i,1) = pp(i,2); % dp'/dy = 0 at y=0 for x < L/4 and x > L/2.
    else
        % dp'/dx = (rho_0*u_0^2*Fxx + dp'/dy)/(dF/dx) at y=0 for L/4 < x < L/2.
        % Fx = (F(i+1) - F(i)) / (dx);
        % Fxx = (F(i+1) - 2*F(i) + F(i-1))/(dx^2);
        % ppy = (pp(i,2) - pp(i,1)) / dy;
        % % TODO: Small Fx could be blowing up the solution!
        % pp(i+1,1) = pp(i-1,1) + (rho_0*u_0^2*Fxx + ppy)/(2*dx*Fx);
        
        % dp'/dy = -rho_0*u_0^2*Fxx
        Fxx = (F(i+1) - 2*F(i) + F(i-1))/(dx^2);
        pp(i,1) = pp(i,2) + rho_0*u_0^2*Fxx*2*dy;
    end
end

% Calculate time derivatives for p' and q' at every interior grid point.
dppdt = zeros(N+1, M+1);
dqpdt = zeros(N+1, M+1);
for i=2:N
    for j=2:M
        dppdt(i,j) = qp(i,j) - u_0 * (pp(i+1,j) - pp(i-1,j))/(2*dx);
        dqpdt(i,j) = (c_0^2) * ((pp(i+1,j) - 2*pp(i,j) + pp(i-1,j)) / (dx^2)) ...
                   + (c_0^2) * ((pp(i,j+1) - 2*pp(i,j) + pp(i,j-1)) / (dy^2)) ...
                   - u_0 * ((qp(i+1,j) - qp(i-1,j)) / (2*dx));
    end
end

dydt = [reshape(dppdt, [(N+1)*(M+1), 1]); reshape(dqpdt, [(N+1)*(M+1), 1])];

end