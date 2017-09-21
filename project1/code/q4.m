function [pp, x, y] = q4(dt, nstep, M, F, ppinitial, qpinitial)
% Defining physical constants.
rho_0 = 1;
c_0 = 1;
u_0 = sqrt(2);
L = 6;
H = 2;
theta = pi/4;

N = 4*(size(F,1)-1); % From the fact that F is of size (N/4 + 1, 1).

% I call them dx and dy but I really mean dx^tilde and dy^tilde. I just
% don't want to keep writing dx_tilde and dy_tilde.
dx = L/N;
dy = H/M;

% Convert F into a vector of size (N+1, 1) so I can index it with i along
% with, e.g. pp(i,j).
F = [zeros(N/4, 1); F; zeros(N/2, 1)];

% Initializing p' and q'.
pp = ppinitial;
qp = qpinitial;

% Using ode45.
y0 = [reshape(pp, [(N+1)*(M+1), 1]); reshape(qp, [(N+1)*(M+1), 1])];
[t, y] = ode45(@(t,y) q4odefun(t, y, N, M, dx, dy, F), [0 nstep*dt], y0, odeset('Stats', 'on'));
pp = reshape(y(1000, 1:(N+1)*(M+1)), [N+1, M+1]);
qp = reshape(y(1000, (N+1)*(M+1)+1:end), [N+1, M+1]);

fprintf('t=%f s to t=%f s\n', t(1), t(end));

x = repmat(linspace(0, L, N+1), M+1, 1);
% y = repmat(linspace(0, H, M+1)', 1, N+1);

yy = repmat(linspace(0, H, M+1)', 1, N+1);
for i = 1:10:size(y,1)
    ppp = reshape(y(i, 1:(N+1)*(M+1)), [N+1, M+1]);
    surf(x,yy,ppp');
    shading interp;
    view(2);
    colorbar;
    % contourf(x,yy,ppp');
    drawnow;
end

y = repmat(linspace(0, H, M+1)', 1, N+1);
end