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

end