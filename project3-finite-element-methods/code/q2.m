function [M2b, ux_N, uy_N] = q2()
%% Global variables
global E nu alpha beta dg_ijdx dg_ijdy dg_ijdx_gen dg_ijdy_gen

%% Physical problem parameters
E = 1.18e11; % Young's modulus for titanium Ti-6Al-2Nb-1Ta-0.8Mo [N/m^2]
nu = 0.31;   % Poisson's ratio for Titanium Ti-6Al-2Nb-1Ta-0.8Mo
alpha = E/(1-nu^2);
beta = (1-nu)/2;
% syms E nu alpha beta; % uncomment to calculate M analytically

%% Function handles for the bilinear basis functions and their derivatives
% Cell array with function handles for the four g_ij basis functions.
% i,j = 1 corresponds to the - subscript. i,j = 2 to the + subscript.
% E.g. g_ij{2,1} returns the g_+- function handle.
g_ij = {@(x,y) (1+x).*(1+y)./4, @(x,y) (1+x).*(1-y)./4;  % g_--, g_-+
        @(x,y) (1-x).*(1+y)./4, @(x,y) (1-x).*(1-y)./4}; % g_+-, g_++

% Generalized basis functions over the domain [x0-a,x0+a]x[y0-b,y0+b]
g_ij_gen = {@(x,y,x0,y0,a,b)  (x-(x0-a))*(y-(y0-b))/(4*a*b), ...
            @(x,y,x0,y0,a,b) -(x-(x0+a))*(y-(y0-b))/(4*a*b);
            @(x,y,x0,y0,a,b) -(x-(x0-a))*(y-(y0+b))/(4*a*b), ...
            @(x,y,x0,y0,a,b)  (x-(x0+a))*(y-(y0+b))/(4*a*b)};
    
% Function handles for the x-derivatives of the basis functions g_ij.
% Note: Returned function handles are still multivariate (x,y) to remain
% compatible with our 2D Gauss-Legendre quadrature function, which
% integrates them over the unit square.
dg_ijdx = {@(x,y)  (1+y)/4, @(x,y)  (1-y)/4;  % dg_--/dx, dg_-+/dx
           @(x,y) -(1+y)/4, @(x,y) -(1-y)/4}; % dg_+-/dx, dg_++/dx

% x derivative of the generalized basis functions
dg_ijdx_gen = {@(x,y,x0,y0,a,b)  (y-(y0-b))/(4*a*b), ...
               @(x,y,x0,y0,a,b) -(y-(y0-b))/(4*a*b);
               @(x,y,x0,y0,a,b) -(y-(y0+b))/(4*a*b), ...
               @(x,y,x0,y0,a,b)  (y-(y0+b))/(4*a*b)};

% Function handles for the y-derivatives of the basis functions g_ij.
dg_ijdy = {@(x,y) (1+x)/4, @(x,y) -(1+x)/4;  % dg_--/dy, dg_-+/dy
           @(x,y) (1-x)/4, @(x,y) -(1-x)/4}; % dg_+-/dy, dg_++/dy

% y derivatives of the generalized basis functions
dg_ijdy_gen = {@(x,y,x0,y0,a,b)  (x-(x0-a))/(4*a*b), ...
               @(x,y,x0,y0,a,b) -(x-(x0+a))/(4*a*b);
               @(x,y,x0,y0,a,b) -(x-(x0-a))/(4*a*b), ...
               @(x,y,x0,y0,a,b)  (x-(x0+a))/(4*a*b)};
       
%% Problem 2(a): single element (unit square) beam under stress
% Element parameters
x0 = 0;
y0 = 0;
a = 1;
b = 1;

% Build M matrix block by block.
M = zeros(8,8);
% M = sym('M', [8 8]);  % uncomment to calculate M analytically
for gamma=1:4
    for delta=1:4
        M(2*gamma-1:2*gamma, 2*delta-1:2*delta) = B(gamma, delta);
    end
end

% Build right hand side vector containing the forces.
F_delta = 8e4; % [N/m^2]

b = zeros(8,1);
g_11 = g_ij{1,1};
b(2) = -F_delta * g_11(1,0);
g_12 = g_ij{1,2};
b(4) = -F_delta * g_12(1,0);
g_21 = g_ij{2,1};
b(6) = -F_delta * g_21(1,0);
g_22 = g_ij{2,2};
b(8) = -F_delta * g_22(1,0);

% Imposing that u_x = 0 and u_y = 0 for x=-1, -1<=y<=1 is equivalent to imposing
% that a_+-, a_++, b+-, and b_++ all equal zero. So we are just interested in
% solving the top left 4x4 system of linear equations.
M2a = M(1:4,1:4);
b2a = b(1:4);
a2a = M2a \ b2a;

% Calculate and plot deflection fields from the a_ij, b_ij coefficients.
ux = @(x,y) a2a(1).*g_ij{1,1}(x,y) + a2a(3).*g_ij{1,2}(x,y);
uy = @(x,y) a2a(2).*g_ij{1,1}(x,y) + a2a(4).*g_ij{1,2}(x,y);
x = repmat(linspace(-1,1,100),[100,1]);
y = repmat(linspace(-1,1,100),[100,1]);
ux_xy = ux(x,y);
uy_xy = uy(x,y);

figure;
surf(x, y', ux_xy); shading interp; colorbar; title('u_x(x,y)');

figure;
surf(x, y', uy_xy); shading interp; colorbar; title('u_y(x,y)');

figure;
plot(x(end/2,:), ux_xy(end/2,:)', x(end/2,:), uy_xy(end/2,:)');
legend('u_x(x,0)', 'u_y(x,0)');

%% Problem 2(b): rectangular beam (composed of N elements) under stress
N = 4;

% Assemble the block-diagonal matrix.
M2b = zeros(4*(N+1), 4*(N+1));
for n=2:N
    x0 = 2*n-1;
    y0 = 0;
    a = 1;
    b = 1;
    
    M2 = zeros(8,8);
    for gamma=1:4
        for delta=1:4
            M2(2*gamma-1:2*gamma, 2*delta-1:2*delta) = B_gen(gamma, delta, x0, y0, a, b);
        end
    end
    
    M2b(4*n-3:4*n+4, 4*n-3:4*n+4) = fliplr(flipud(M2));
end

% Imposing the u_x = u_y = 0 BC at the left wall.
M2b(1:4,1:4) = zeros(4,4);
M2b(5:8,1:4) = zeros(4,4);
M2b(1:4,5:8) = zeros(4,4);

b2b = zeros(4*N+4,1);
b2b(end) = -F_delta * g_11(2*N,0);
b2b(end-2) = -F_delta * g_12(2*N,0);
b2b(end-4) = -F_delta * g_21(2*N,0);
b2b(end-6) = -F_delta * g_22(2*N,0);

% b2b(end-6) = -F_delta * g_11(2*N,0);
% b2b(end-4) = -F_delta * g_12(2*N,0);
% b2b(end-2) = -F_delta * g_21(2*N,0);
% b2b(end) = -F_delta * g_22(2*N,0);

a2b = zeros(4*N+4,1);
a2b(5:end) = M2b(5:end, 5:end) \ b2b(5:end);
% a2b = M2b \ b2b;

ux_N = [];
uy_N = [];

for n=1:N
    a_n = a2b(4*n-3:4*n+4);
    aij = a_n(1:2:7);
    bij = a_n(2:2:8);
 
    x0 = 2*n-1;
    y0 = 0;
    a = 1;
    b = 1;
    ux = @(x,y) aij(4).*g_ij_gen{1,1}(x,y,x0,y0,a,b) ...
              + aij(3).*g_ij_gen{1,2}(x,y,x0,y0,a,b) ...
              + aij(2).*g_ij_gen{2,1}(x,y,x0,y0,a,b) ...
              + aij(1).*g_ij_gen{2,2}(x,y,x0,y0,a,b);
    uy = @(x,y) bij(4).*g_ij_gen{1,1}(x,y,x0,y0,a,b) ...
              + bij(3).*g_ij_gen{1,2}(x,y,x0,y0,a,b) ...
              + bij(2).*g_ij_gen{2,1}(x,y,x0,y0,a,b) ...
              + bij(1).*g_ij_gen{2,2}(x,y,x0,y0,a,b);
    
    x = repmat(linspace(2*n,2*n+2,100), [100,1]);
    y = repmat(linspace(-1,1,100), [100,1]);
    
    ux_N = [ux_N; ux(x,y)'];
    uy_N = [uy_N; uy(x,y)'];
end

x = linspace(0,2*N,100*N);
y = linspace(-1,1,100);
figure; surf(x,y,ux_N'); title('2(b) u_x'); shading interp; colorbar;
figure; surf(x,y,uy_N'); title('2(b) u_y'); shading interp; colorbar;
end

%% Gauss-Legendre quadrature integration
function integral = gauss_legendre_quadrature_2D(f, x1, x2, y1, y2)
% Integrate the function f(x,y) over the domain [x1,x2]x[y1,y2] using
% Gaussian quadrature with a unit weighing function \omega(x,y)=1 so the
% associated polynomials are the Legendre polynomials and thus the method
% is known as Gauss-Legendre quadrature with two quadrature points for each
% spatial dimension as we only need to integrate bilinear functions.
    scaled_f = @(x,y) f((x2-x1)*x/2 + (x1+x2)/2, (y2-y1)*y/2 + (y1+y2)/2);
    integral = (x2-x1)*(y2-y1)/4 * ...
               ( scaled_f(-1/sqrt(3), -1/sqrt(3)) ...
               + scaled_f(-1/sqrt(3), 1/sqrt(3)) ...
               + scaled_f(1/sqrt(3), -1/sqrt(3)) ...
               + f(1/sqrt(3), 1/sqrt(3)) );
end

%% Interface for looking up derivatives of the bilinear basis functions.
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

function dg = dg_gen(i,j,m)
% Returns the function handle corresponding to the m derivative of the g_ij
% basis function, that is dg_ij/dx_m. i,j can take on the values {1,2} where 1
% corresponds to the - subscript, 2 to the + subscript. m can take on the string
% values {'x','y'}.

    assert(i==1 || i==2, 'Invalid value for i, not 1 (-) or 2 (+)!')
    assert(j==1 || j==2, 'Invalid value for j, not 1 (-) or 2 (+)!')
    assert(m=='x' || m=='y', 'Invalid value for m, not ''x'' or ''y''!')

    global dg_ijdx_gen dg_ijdy_gen

    if m == 'x'
        dg = dg_ijdx_gen{i,j};
    elseif m == 'y'
        dg = dg_ijdy_gen{i,j};
    end
end

%% Implementation of the rank-6 G_ijkl^mn tensor-like symbol.
function G = G(i,j,k,l,m,n)
    global alpha
    assert(i==1 || i==2, 'Invalid value for i, not 1 (-) or 2 (+)!')
    assert(j==1 || j==2, 'Invalid value for j, not 1 (-) or 2 (+)!')
    assert(k==1 || k==2, 'Invalid value for k, not 1 (-) or 2 (+)!')
    assert(l==1 || l==2, 'Invalid value for l, not 1 (-) or 2 (+)!')
    assert(m=='x' || m=='y', 'Invalid value for m, not ''x'' or ''y''!')
    assert(n=='x' || n=='y', 'Invalid value for m, not ''x'' or ''y''!')

    dg1 = dg(i,j,m);
    dg2 = dg(k,l,n);
    integrand = @(x,y) alpha*dg1(x,y)*dg2(x,y);
    G = gauss_legendre_quadrature_2D(integrand, -1, 1, -1, 1);
end

function G = G_gen(i,j,k,l,m,n,x0,y0,a,b)
    global alpha
    assert(i==1 || i==2, 'Invalid value for i, not 1 (-) or 2 (+)!')
    assert(j==1 || j==2, 'Invalid value for j, not 1 (-) or 2 (+)!')
    assert(k==1 || k==2, 'Invalid value for k, not 1 (-) or 2 (+)!')
    assert(l==1 || l==2, 'Invalid value for l, not 1 (-) or 2 (+)!')
    assert(m=='x' || m=='y', 'Invalid value for m, not ''x'' or ''y''!')
    assert(n=='x' || n=='y', 'Invalid value for m, not ''x'' or ''y''!')

    dg1 = dg_gen(i,j,m);
    dg2 = dg_gen(k,l,n);
    integrand = @(x,y) alpha*dg1(x,y,x0,y0,a,b)*dg2(x,y,x0,y0,a,b);
    G = gauss_legendre_quadrature_2D(integrand, x0-a, x0+a, y0-b, y0+b);
end

%% Generates and returns the 2x2 block matrices B_gamma,delta.
function B = B(gamma, delta)
    global nu beta

    four2two = {[1,1], [1,2], [2,1], [2,2]};
    [i,j] = deal(four2two{gamma}(1), four2two{gamma}(2));
    [k,l] = deal(four2two{delta}(1), four2two{delta}(2));

    B_11 = G(i,j,k,l,'x','x') + beta*G(i,j,k,l,'y','y');
    B_12 = nu*G(i,j,k,l,'x','y') + beta*G(i,j,k,l,'y','x');
    B_21 = beta*G(i,j,k,l,'x','y') + nu*G(i,j,k,l,'y','x');
    B_22 = beta*G(i,j,k,l,'x','x') + G(i,j,k,l,'y','y');

    B = [B_11, B_12; B_21, B_22];
end

function B = B_gen(gamma, delta,x0,y0,a,b)
    global nu beta

    four2two = {[1,1], [1,2], [2,1], [2,2]};
    [i,j] = deal(four2two{gamma}(1), four2two{gamma}(2));
    [k,l] = deal(four2two{delta}(1), four2two{delta}(2));
    
    G2 = @(i,j,k,l,m,n) G_gen(i,j,k,l,m,n,x0,y0,a,b);

    B_11 = G2(i,j,k,l,'x','x') + beta*G2(i,j,k,l,'y','y');
    B_12 = nu*G2(i,j,k,l,'x','y') + beta*G2(i,j,k,l,'y','x');
    B_21 = beta*G2(i,j,k,l,'x','y') + nu*G2(i,j,k,l,'y','x');
    B_22 = beta*G2(i,j,k,l,'x','x') + G2(i,j,k,l,'y','y');

    B = [B_11, B_12; B_21, B_22];
end