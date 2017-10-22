function q2()

    % Integrate the function f(x,y) over the unit square [-1,1]x[-1,1]
    % using Gaussian quadrature with a unit weighing function \omega(x,y)=1
    % so the associated polynomials are the Legendre polynomials and thus
    % the method is known as Gauss-Legendre quadrature with two evaluation
    % points for each dimension.
    function integral = gauss_legendre_quadrature_2D(f)
        integral = f(-1/sqrt(3), -1/sqrt(3)) + f(-1/sqrt(3), 1/sqrt(3)) ...
                   + f(1/sqrt(3), -1/sqrt(3)) + f(1/sqrt(3), 1/sqrt(3));
    end

g_ij = {@(x,y) (1+x)*(1+y)/4, @(x,y) (1+x)*(1-y)/4;  % g_--, g_-+
        @(x,y) (1-x)*(1+y)/4, @(x,y) (1-x)*(1-y)/4}; % g_+-, g_++

dg_ijdx = {@(x,y) (1+y)/4,  @(x,y) (1-y)/4;   % dg_--/dx, dg_-+/dx
           @(x,y) -(1+y)/4, @(x,y) -(1-y)/4}; % dg_+-/dx, dg_++/dx

dg_ijdy = {@(x,y) (1+x)/4, @(x,y) -(1+x)/4;  % dg_--/dy, dg_-+/dy
           @(x,y) (1-x)/4, @(x,y) -(1-x)/4}; % dg_+-/dy, dg_++/dy

end