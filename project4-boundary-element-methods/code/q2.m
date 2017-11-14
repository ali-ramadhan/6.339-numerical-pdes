n = 15;
theta = linspace(-pi,pi,n+1);
w = 2*pi/n;
A = -pi*eye(n+1) - ones(n+1,n+1)*(w/2);
b = 1 ./ (3 + 2*cos(theta') + cos(2*theta'));
sigma = A \ b;