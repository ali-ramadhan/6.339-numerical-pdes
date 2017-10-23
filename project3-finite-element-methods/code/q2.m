function q2()
global E nu dg_ijdx dg_ijdy

% Physical parameters
E = 1.18e11; % Young's modulus for titanium Ti-6Al-2Nb-1Ta-0.8Mo [N/m^2]
nu = 0.31;   % Poisson's ratio for Titanium Ti-6Al-2Nb-1Ta-0.8Mo

% Cell array with function handles for the four g_ij basis functions.
% i,j = 1 corresponds to the - subscript. i,j = 2 to the + subscript.
% E.g. g_ij{2,1} returns the g_+- function handle.
g_ij = {@(x,y) (1+x).*(1+y)./4, @(x,y) (1+x).*(1-y)./4;  % g_--, g_-+
        @(x,y) (1-x).*(1+y)./4, @(x,y) (1-x).*(1-y)./4}; % g_+-, g_++

% Function handles for the x-derivatives of the basis functions g_ij.
% Note: Returned function handles are still multivariate (x,y) to remain
% compatible with our 2D Gauss-Legendre quadrature function, which
% integrates them over the unit square.
dg_ijdx = {@(x,y)  (1+y)/4, @(x,y)  (1-y)/4;  % dg_--/dx, dg_-+/dx
           @(x,y) -(1+y)/4, @(x,y) -(1-y)/4}; % dg_+-/dx, dg_++/dx

% Function handles for the y-derivatives of the basis functions g_ij.
dg_ijdy = {@(x,y) (1+x)/4, @(x,y) -(1+x)/4;  % dg_--/dy, dg_-+/dy
           @(x,y) (1-x)/4, @(x,y) -(1-x)/4}; % dg_+-/dy, dg_++/dy

M = zeros(8,8);
% M = [B(1,1), B(1,2), B(1,3), B(1,4);
%      B(2,1), B(2,2), B(2,3), B(2,4);
%      B(3,1), B(3,2), B(3,3), B(3,4);
%      B(4,1), B(4,2), B(4,3), B(4,4)];
for gamma=1:4
    for delta=1:4
        M(2*gamma-1:2*gamma, 2*delta-1:2*delta) = B(gamma, delta);
    end
end

M
det(M)
cond(M)

F_delta = 8e4; % [N/m^2]
alpha = E/(1-nu^2);

b = zeros(8,1);
g_11 = g_ij{1,1};
b(2) = -(F_delta/alpha) * g_11(1,0);
g_12 = g_ij{1,2};
b(4) = -(F_delta/alpha) * g_12(1,0);
g_21 = g_ij{2,1};
b(6) = -(F_delta/alpha) * g_21(1,0);
g_22 = g_ij{2,2};
b(8) = -(F_delta/alpha) * g_22(1,0);

b

a = M \ b

ux = @(x,y) a(1).*g_ij{1,1}(x,y) + a(3).*g_ij{1,2}(x,y) + a(5).*g_ij{2,1}(x,y) + a(7).*g_ij{2,2}(x,y);
uy = @(x,y) a(2).*g_ij{1,1}(x,y) + a(4).*g_ij{1,2}(x,y) + a(6).*g_ij{2,1}(x,y) + a(8).*g_ij{2,2}(x,y);
x = repmat(linspace(-1,1,10),[10,1]);
y = repmat(linspace(-1,1,10),[10,1]);
ux_xy = ux(x,y);
uy_xy = uy(x,y);
% surf(x,y',ux_xy);
% heatmap(ux_xy);
% figure;
% heatmap(uy_xy);
% figure;
% ux_xy(1,:)
% uy_xy(1,:)
plot(x(end,:), ux_xy(end,:)', x(end,:), uy_xy(end,:)');

end

function integral = gauss_legendre_quadrature_2D(f)
% Integrate the function f(x,y) over the unit square [-1,1]x[-1,1] using
% Gaussian quadrature with a unit weighing function \omega(x,y)=1 so the
% associated polynomials are the Legendre polynomials and thus the method is
% known as Gauss-Legendre quadrature with two evaluation points for each spatial
% dimension.
    
integral = f(-1/sqrt(3), -1/sqrt(3)) + f(-1/sqrt(3), 1/sqrt(3)) ...
           + f(1/sqrt(3), -1/sqrt(3)) + f(1/sqrt(3), 1/sqrt(3));
end

function dg = dg(i,j,m)
% Returns the function handle corresponding to the m derivative of the g_ij
% basis function, that is dg_ij/dx_m. i,j can take on the values {1,2} where 1
% corresponds to the - subscript, 2 to the + subscript. m can take on the string
% values {'x','y'}.

assert(i==1 || i==2, 'Invalid value for i, not 1 (-) or 2 (+)!')
assert(j==1 || j==2, 'Invalid value for j, not 1 (-) or 2 (+)!')
assert(m=='x' || m=='y', 'Invalid value for m, not ''x'' or ''y''!')

global dg_ijdx dg_ijdy

if m == 'x'
    dg = dg_ijdx{i,j};
elseif m == 'y'
    dg = dg_ijdy{i,j};
end
end

function G = G(i,j,k,l,m,n)
assert(i==1 || i==2, 'Invalid value for i, not 1 (-) or 2 (+)!')
assert(j==1 || j==2, 'Invalid value for j, not 1 (-) or 2 (+)!')
assert(k==1 || k==2, 'Invalid value for k, not 1 (-) or 2 (+)!')
assert(l==1 || l==2, 'Invalid value for l, not 1 (-) or 2 (+)!')
assert(m=='x' || m=='y', 'Invalid value for m, not ''x'' or ''y''!')
assert(n=='x' || n=='y', 'Invalid value for m, not ''x'' or ''y''!')

global E nu
alpha = E/(1-nu^2);

dg1 = dg(i,j,m);
dg2 = dg(k,l,n);
integrand = @(x,y) alpha*dg1(x,y)*dg2(x,y);
G = gauss_legendre_quadrature_2D(integrand);
end

function B = B(gamma, delta)
global nu
beta = (1-nu)/2;

four2two = {[1,1], [1,2], [2,1], [2,2]};
[i,j] = deal(four2two{gamma}(1), four2two{gamma}(2));
[k,l] = deal(four2two{delta}(1), four2two{delta}(2));

B_11 = G(i,j,k,l,'x','x') + beta*G(i,j,k,l,'y','y');
B_12 = nu*G(i,j,k,l,'x','y') + beta*G(i,j,k,l,'y','x');
B_21 = beta*G(i,j,k,l,'x','y') + nu*G(i,j,k,l,'y','x');
B_22 = beta*G(i,j,k,l,'x','x') + G(i,j,k,l,'y','y');

B = [B_11, B_12; B_21, B_22];
end