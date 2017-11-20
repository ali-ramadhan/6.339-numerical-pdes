%% 2(a) Solving the exterior Neumann problem using the Nystrom method and plotting sigma
N = [20 50 100 500];

for i=1:size(N,2)
    n = N(i);
    
    theta = linspace(-pi,pi,n+1);
    theta = theta(1:end-1); % Get rid of overlapping endpoint.
    w = 2*pi/n;
    A = -pi*eye(n) - ones(n,n)*(w/2);
    b = 1 ./ (3 + 2*cos(theta') + cos(2*theta'));
    sigma = A \ b;
       
    % Evaluate and plot the solution on a finer grid
    n_eval = 1000;
    theta_eval = linspace(min(theta), max(theta), n_eval);
    sigma_eval = zeros(n_eval,1);
    count = 2;
    for j=1:n_eval
        if theta_eval(j) <= theta(count)
            %count stays same
        else
            count=count+1;
        end
        sigma_eval(j) = sigma(count-1);
    end

    plot(theta_eval, sigma_eval, 'LineWidth', 2, 'DisplayName', strcat('$$n = ', num2str(n), '$$'))
    hold on;
end

% plot exact solution
I = -(1/6) * sqrt(3 + 2*sqrt(3));
f = @(x) 1 ./ (3 + 2.*cos(x) + cos(2.*x));
theta = linspace(-pi, pi, 1000);
sigma_exact = @(x) -(1/pi) * (I/2 + f(x));
sigma_exact_eval = sigma_exact(theta);
plot(theta, sigma_exact_eval, 'LineWidth', 2, 'DisplayName', 'Exact');

legend('show');


%% 2(b) Investigating the convergence of the Nystrom method for the exterior Neumann problem
figure;

N = linspace(1,500,500);
error = zeros(500,1);

for i=1:size(N,2)
    n = N(i);
    
    theta = linspace(-pi,pi,n+1);
    theta = theta(1:end-1); % Get rid of overlapping endpoint.
    w = 2*pi/n;
    A = -pi*eye(n) - ones(n,n)*(w/2);
    b = 1 ./ (3 + 2*cos(theta') + cos(2*theta'));
    sigma = A \ b;
    
    for j=1:n
        error(i) = error(i) + abs(sigma(i) - sigma_exact(theta(i)));
    end
end

loglog(N, error, 'LineWidth', 2);

%% 2(d) Solving the exterior problem on an elliptical boundary
% a = 1.5; % semi-minor axis
% b = 1.2; % semi-major axis
% e = sqrt(1 - (b^2/a^2)); % eccentricity
% p = a*ellipticE(2*pi, e); % perimeter of an ellipse
% 
% n = 10;
% theta = linspace(-pi,pi,10*n);
% theta_eval = zeros(n,1);
% x = zeros(n,1);
% y = zeros(n,1);
% nx = zeros(n,1);
% ny = zeros(n,1);
% bi = zeros(n,1); % RHS vector
% 
% index = 1;
% for i=1:10*n
%     t = atan((a/b) * tan(theta(i) + pi));
%     arc_length = a*ellipticE(t, e);
%     if arc_length >= ((index-1)/n)*p
%         theta_eval(index) = theta(i);
%         x(index) = a*cos(theta(i));
%         y(index) = b*sin(theta(i));
%         nx(index) = -b*sin(theta(i));
%         ny(index) = a*cos(theta(i));
%         bi(index) = 1 / (3 + 2*cos(theta(i)) + cos(2*theta(i)));
%         index = index + 1;
%     end
% end
% 
% theta_eval
% 
% w = p/n;
% A = zeros(n);
% 
% for i=1:n
%     for j=1:n
%         if i==j
%             A(i,j) = -pi;
%         else
%             up = (x(i) - x(j))*nx(i) + (y(i) - y(j))*ny(i);
%             down = (x(i) - x(j))^2 + (y(i) - y(j))^2;
%             A(i,j) = -w * (up/down);
%         end
%     end
% end
% 
% A
% 
% figure;
% imagesc(A);
% colorbar();
% 
% sigma = A \ bi;
% 
% figure;
% plot(theta_eval, sigma);