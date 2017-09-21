function dydt = q4odefun(t, y, N, M, dx, dy, F)
% Redefining problem parameters here to avoid a massive argument list.
% Note: all constants are unitless, or rather, non-dimensional.
rho_0 = 1;
c_0 = 1;
u_0 = sqrt(2);
L = 6;
H = 2;
theta = pi/4;

% Reshape y vector into p',q' matrices.
pp = reshape(y(1:(N+1)*(M+1)), N+1, M+1);
qp = reshape(y((N+1)*(M+1)+1:end), N+1, M+1);

% Imposing boundary conditions at each wall.
pp(1,:) = zeros(1,M+1); % p=p_0 (p'=0) at x=0
pp(N+1,:) = pp(N,:);    % dp'/dx = 0 at x=L
qp(1,:) = zeros(1,M+1); % q'=0 at x=0
qp(N+1,:) = qp(N,:);    % dq'/dx 0 at x=L

% Solve system of linear equations for the y=0 and y=H boundary conditions.
A0 = - (3*sec(theta) / (2*dy)) * diag(ones(N+1,1)) ...
    + (tan(theta) / (2*dx)) * diag(ones(N,1), -1) ...
    - (tan(theta) / (2*dx)) * diag(ones(N,1), 1);
Fxx = [0; diff(F,2)/dx^2; 0];
b0 = - sec(theta)*(4*pp(:,2) - pp(:,3)) / (2*dy) - rho_0*(u_0)^2*Fxx;
pp(:,1) = A0\b0;

AH = (3*sec(theta) / (2*dy)) * diag(ones(N+1,1)) ...
    + (tan(theta) / (2*dx)) * diag(ones(N,1), -1) ...
    - (tan(theta) / (2*dx)) * diag(ones(N,1), 1);
bH = sec(theta)*(4*pp(:,M) - pp(:,M-1)) / (2*dy);
pp(:,M+1) = AH\bH;

% Calculate time derivatives for p' and q' at every interior grid point.
dppdt = zeros(N+1, M+1);
dqpdt = zeros(N+1, M+1);
for i=2:N
    for j=2:M
        dppdt(i,j) = qp(i,j) - u_0 * (pp(i+1,j) - pp(i-1,j))/(2*dx);
        
        dqpdx = (qp(i+1,j) - qp(i-1,j)) / (2*dx);
        dppdxx = (pp(i+1,j) - 2*pp(i,j) + pp(i-1,j)) / (dx^2);
        dppdyy = (pp(i,j+1) - 2*pp(i,j) + pp(i,j-1)) / (dy^2);
        dppdxy = (pp(i+1,j+1) - pp(i+1,j-1) - pp(i-1,j+1) ...
            + pp(i-1,j-1)) / (4*dx*dy);
        
        dqpdt(i,j) = (c_0^2)*(sec(theta)^2)*(dppdxx + dppdyy) ...
                   - (c_0^2)*sec(theta)*tan(theta)*dppdxy ...
                   - u_0*dqpdx;
    end
end

dydt = [reshape(dppdt, [(N+1)*(M+1), 1]); reshape(dqpdt, [(N+1)*(M+1), 1])];
end