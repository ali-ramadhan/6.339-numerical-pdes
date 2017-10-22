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

    function G = G(i,j,k,l,m,n)
        % Returns a function corresponding to the m derivative of the g_ij
        % function, that is dg_ij/dx_m, where i,j = 1,2 where {1,2}->{-,+}
        % and m = x,y.
        function dg = dg(i,j,m)
            if m=='x'
                dg = dg_ijdx{i,j};
            elseif m=='y'
                dg = dg_ijdy{i,j};
            end 
        end
        
        E = 1;
        nu = 0;
        alpha = E/(1-nu^2);
        dg1 = dg(i,j,m);
        dg2 = dg(k,l,n);
        integrand = @(x,y) alpha*dg1(x,y)*dg2(x,y);
        G = gauss_legendre_quadrature_2D(integrand);
    end
       
end