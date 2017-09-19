function [pp, x, y] = q2(dt, nstep, M, F, ppinitial, qpinitial)

global N dx dy

rho_0 = 1;
c_0 = 1;
u_0 = sqrt(2);
L = 6;
H = 2;

N = 4*(size(F,1)-1); % From the fact that F is of size (N/4 + 1, 1).

dx = L/N;
dy = H/M;

% Convert F into a vector of size (N+1, 1) so I can index it with i along
% with, e.g. pp(i,j).
F = [zeros(N/4, 1); F; zeros(N/2, 1)];
% fprintf('size(F,1) = %d\n', size(F,1));
% fprintf('N = %d\n', N);

% Initializing p' and q'.
pp = ppinitial;
qp = qpinitial;

% Using ode45.
y0 = [reshape(pp, [(N+1)*(M+1), 1]); reshape(qp, [(N+1)*(M+1), 1])];
[t, y] = ode45(@(t,y) q2odefun(t, y, M, F), [0 nstep*dt], y0, odeset('Stats', 'on'));
pp = reshape(y(end, 1:(N+1)*(M+1)), [N+1, M+1]);
qp = reshape(y(end, (N+1)*(M+1)+1:end), [N+1, M+1]);

fprintf('t=%f s to t=%f s\n', t(1), t(end));

% Euler's method
% for n=1:nstep
%     % Imposing boundary conditions.    
%     for j=1:M+1
%         pp(1,j) = 0; % p=p_0 (p'=0) at x=0
%         pp(N+1,j) = pp(N,j); % dp'/dx = 0 at x=L
%         
%         qp(1,j) = 0; % q'=0 at x=0
%         qp(N+1,j) = qp(N,j); % dq'/dx 0 at x=L
%     end
%     
%     for i=1:N+1
%         pp(i,M+1) = pp(i,M); % dp'/dy = 0 at y=H
%         
%         fprintf('i = %d\n', i);
%         fprintf('x = %d\n', i*dx);
%         if i*dx < L/4 || i*dx > L/2
%             pp(i,1) = pp(i,2); % dp'/dy = 0 at y=0 for x < L/4 and x > L/2.
%         else
%             % dp'/dx = (rho_0*u_0^2*Fxx + dp'/dy)/(dF/dx) at y=0 for L/4 < x < L/2.
%             Fx = (F(i+1) - F(i)) / (dx);
%             Fxx = (F(i+1) - 2*F(i) + F(i-1))/(dx^2);
%             ppy = (pp(i,2) - pp(i,1)) / dy;
%             pp(i+1,1) = pp(i-1,1) + (rho_0*u_0^2*Fxx + ppy)/(2*dx*Fx);
%         end
%     end
%     
%     % Calculate time derivatives for p' and q' at every grid point.
%     % We will not worry about calculating derivatives at the boundaries.
%     for i=2:N
%         for j=2:M
%             dppdt(i,j) = qp(i,j) - u_0 * (pp(i+1,j) - pp(i-1,j))/(2*dx);
%             dqpdt(i,j) = (c_0^2) * ((pp(i+1,j) - 2*pp(i,j) + pp(i-1,j)) / (dx^2)) ...
%                        + (c_0^2) * ((pp(i,j+1) - 2*pp(i,j) + pp(i,j-1)) / (dy^2)) ...
%                        - u_0 * ((qp(i+1,j) - qp(i-1,j)) / (2*dx));
%         end
%     end
% 
%     % Calculate new values of p' and q' at every grid point (Euler's method).
%     % The boundary is taken care of by the boudary conditions.
%     for i=2:N
%         for j=2:M
%             pp(i,j) = pp(i,j) + dppdt(i,j)*dt;
%             qp(i,j) = qp(i,j) + dqpdt(i,j)*dt;
%         end
%     end
%     
% end

x = repmat(linspace(0, L, N+1), M+1, 1);
yy = repmat(linspace(0, H, M+1)', 1, N+1);

for i = 1:1:size(y,1)
    ppp = reshape(y(i, 1:(N+1)*(M+1)), [N+1, M+1]);
    %surf(x,yy,ppp');
    %shading interp;
    %view(2);
    %colorbar;
    contour(x,yy,ppp');
    drawnow;
end

y = repmat(linspace(0, H, M+1)', 1, N+1);

end