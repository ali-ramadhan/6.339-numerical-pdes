%% Solving the exterior Neumann problem using the Nystrom method and plotting sigma
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
    theta_eval = linspace(min(theta), max(theta), nEval);
    sigma_eval = zeros(nEval,1);
    count = 2;
    for j=1:n_eval
        if theta_eval(j) <= theta(count)
            %count stays same
        else
            count=count+1;
        end
        sigma_eval(j) = sigma(count-1);
    end

    plot(xEval, sigma_eval, 'LineWidth', 2, 'DisplayName', strcat('$$n = ', num2str(n), '$$'))
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


%% Investigating the convergence of the Nystrom method for the exterior Neumann problem
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