function dydt = q2odefun(t, y, N, M, dx, dy, F)
global pp qp

% Redefining problem parameters here to avoid a massive argument list.
% Note: all constants are unitless, or rather, non-dimensional.
rho_0 = 1;
c_0 = 1;
u_0 = sqrt(2);
L = 6;
H = 2;

% Reshape y vector into p' and q' matrices for more readable code.
pp = reshape(y(1:(N+1)*(M+1)), N+1, M+1);
qp = reshape(y((N+1)*(M+1)+1:end), N+1, M+1);

% Imposing boundary conditions at each wall.
pp(1,:) = zeros(1,M+1); % p=p_0 (p'=0) at x=0
pp(N+1,:) = pp(N,:);    % dp'/dx = 0 at x=L
qp(1,:) = zeros(1,M+1); % q'=0 at x=0
qp(N+1,:) = qp(N,:);    % dq'/dx 0 at x=L

pp(:,M+1) = pp(:,M); % dp'/dy = 0 at y=H

for i=1:N+1
    if i*dx < L/4 || i*dx > L/2
        % Imposing dp'/dy = 0 at y=0 for x < L/4 and x > L/2.
        pp(i,1) = pp(i,2);
    else     
        % Imposing dp'/dy = -rho_0*u_0^2*Fxx at y=0 for L/4 <= x <= L/2.
        Fxx = (F(i+1) - 2*F(i) + F(i-1))/(dx^2);

        % Multiply by 2*dy to continue using central difference with
        % second-order truncation error (but now essentially with a ghost
        % point).
        pp(i,1) = pp(i,2) + rho_0*u_0^2*Fxx*dy;
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

dydt = [reshape(dppdt, (N+1)*(M+1), 1); reshape(dqpdt, (N+1)*(M+1), 1)];
end