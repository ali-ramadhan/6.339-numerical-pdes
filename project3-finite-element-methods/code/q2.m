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

g_ij = {@(x,y) (1+x).*(1+y)./4, @(x,y) (1+x).*(1-y)./4;  % g_--, g_-+
        @(x,y) (1-x).*(1+y)./4, @(x,y) (1-x).*(1-y)./4}; % g_+-, g_++

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
        
        E = 1.18e11;  % Young's modulus for Ti-6Al-2Nb-1Ta-0.8Ma [N/m^2]
        nu = 0.31;
        alpha = E/(1-nu^2);
        dg1 = dg(i,j,m);
        dg2 = dg(k,l,n);
        integrand = @(x,y) alpha*dg1(x,y)*dg2(x,y);
        G = gauss_legendre_quadrature_2D(integrand);
    end

four2two = {[1,1], [1,2], [2,1], [2,2]};

    function B = B(gamma, delta)
        function [i,j] = four2two(mu)
            if mu == 1
                i=1; j=1;
            elseif mu == 2
                i=1; j=2;
            elseif mu == 3
                i=2; j=1;
            elseif mu == 4
                i=2; j=2;
            end
        end
        
        nu = 0.31;
        beta = (1-nu)/2;
        [i,j] = four2two(gamma);
        [k,l] = four2two(delta);
        B_11 = G(i,j,k,l,'x','x') + beta*G(i,j,k,l,'y','y');
        B_12 = nu*G(i,j,k,l,'x','y') + beta*G(i,j,k,l,'y','x');
        B_21 = beta*G(i,j,k,l,'x','y') + nu*G(i,j,k,l,'y','x');
        B_22 = beta*G(i,j,k,l,'x','x') + G(i,j,k,l,'y','y');
        B = [B_11, B_12; B_21, B_22];
    end

M = zeros(8,8);
for gamma=1:4
    for delta=1:4
        M(2*gamma-1:2*gamma, 2*delta-1:2*delta) = B(gamma, delta);
    end
end

E = 1.18e11;  % Young's modulus for Ti-6Al-2Nb-1Ta-0.8Ma [N/m^2]
nu = 0.31;
alpha = E/(1-nu^2);

b = zeros(8,1);
% k = 0;
% for i=1:2
%     for j =1:2
%         k = k+1;
%         g = g_ij{i,j};
%         integrand_x = @(x,y) (1/alpha)*g(x,y)*Fx;
%         integrand_y = @(x,y) (1/alpha)*g(x,y)*Fy;
%         b(k) = gauss_legendre_quadrature_2D(integrand_x);
%         b(4+k) = gauss_legendre_quadrature_2D(integrand_y);
%     end
% end

% heatmap(M)

Fx = 0;
Fy = -8e4;

g = g_ij{2,1};
integrand_y = @(x,y) (1/alpha)*g(x,y)*Fy;
b(6) = gauss_legendre_quadrature_2D(integrand_y);
a = M \ b;

ux = @(x,y) a(1).*g_ij{1,1}(x,y) + a(3).*g_ij{1,2}(x,y) + a(5).*g_ij{2,1}(x,y) + a(7).*g_ij{2,2}(x,y);
uy = @(x,y) a(2).*g_ij{1,1}(x,y) + a(4).*g_ij{1,2}(x,y) + a(6).*g_ij{2,1}(x,y) + a(8).*g_ij{2,2}(x,y);
x = repmat(linspace(-1,1,10),[10,1]);
y = repmat(linspace(-1,1,10),[10,1]);
ux_xy = ux(x,y);
uy_xy = uy(x,y);
heatmap(ux_xy);
figure;
heatmap(uy_xy);
figure;
ux_xy(1,:)
uy_xy(1,:)
plot(x(end,:), ux_xy(end,:)', x(end,:), uy_xy(end,:)');

end