function [pp, x, y] = q2(dt, nstep, M, F, ppinitial, qpinitial)
global pp qp

% Defining problem parameters. Note: all constants are unitless, or
% rather, non-dimensional.
rho_0 = 1;
c_0 = 1;
u_0 = sqrt(2);
L = 6;
H = 2;

% From the fact that F is of size (N/4 + 1, 1).
N = 4*(size(F,1)-1);

% Determine step sizes.
dx = L/N;
dy = H/M;

% Convert F into a vector of size (N+1, 1) so I can index it with i along
% with, e.g. pp(i,j).
F = [zeros(N/4, 1); F; zeros(N/2, 1)];

% Initializing p' and q'.
pp = ppinitial;
qp = qpinitial;

% Use ode45 to integrate 
y0 = [reshape(pp, (N+1)*(M+1), 1); reshape(qp, (N+1)*(M+1), 1)];
[t, y] = ode45(@(t,y) q2odefun(t, y, N, M, dx, dy, F), [0 nstep*dt], y0, odeset('Stats', 'on'));

x = repmat(linspace(0, L, N+1), M+1, 1);
y = repmat(linspace(0, H, M+1)', 1, N+1);
end