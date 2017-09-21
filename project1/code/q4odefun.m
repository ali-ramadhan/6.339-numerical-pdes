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
for j=1:M+1
    pp(1,j) = 0; % p=p_0 (p'=0) at x=0
    pp(N+1,j) = pp(N,j); % dp'/dx = 0 at x=L

    qp(1,j) = 0; % q'=0 at x=0
    qp(N+1,j) = qp(N,j); % dq'/dx 0 at x=L
end

% for i=3:N+1
%     % TODO: Why start at i=2?
%     % up = tan(theta)*pp(i-1,M+1)*dy - sec(theta)*pp(i,M)*dx;
%     % down = tan(theta)*dy - sec(theta)*dx;
%     % pp(i,M+1) = up/down;
%     
%     if i*dx < L/4 || i*dx > L/2
%         up = tan(theta)*pp(i-2,1)*dy - sec(theta)*pp(i,3)*dx;
%         down = tan(theta)*dy - sec(theta)*dx;
%         pp(i,1) = up/down;
%     else
%         Fxx = (F(i+1) - 2*F(i) + F(i-1))/(dx^2);
%         up = rho_0*u_0^2*Fxx*dx*dy + tan(theta)*pp(i-1,1)*dy ...
%             + sec(theta)*pp(i+2,2)*dx;
%         down = tan(theta)*dy + sec(theta)*dx;
%         pp(i,1) = up/down;
%     end
% end

% Solve system of linear equations for the y=0 and y=H boundaries.
A = - (3*sec(theta) / (2*dy)) * diag(ones(N+1,1)) ...
    + (tan(theta) / (2*dx)) * diag(ones(N,1), -1) ...
    - (tan(theta) / (2*dx)) * diag(ones(N,1), 1);
bH = - sec(theta)*(4*pp(:,2) - pp(:,3)) / (2*dy);

Fxx = [0; diff(F,2)/dx^2; 0];
b0 = bH - rho_0*(u_0)^2*Fxx;

% ppyH = A\bH;
% pp(2:N,M+1) = ppyH(2:N);
% ppy0 = A\b0;
% pp(2:N,1) = ppy0(2:N);
pp(:,M+1) = A\bH;
pp(:,1) = A\b0;
% pp(1,1) = 0; pp(N+1,1) = pp(N,1);
% pp(1,M+1) = 0; pp(N+1,M+1) = pp(N,M+1);

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