A = [2, 0, 0, 0];
B = [10, pi, 1, 2*pi];
F = {@(x) 3*x, @(x) sin(x), @(x) exp(2*cos(2*pi*x)), @(x) abs(cos(x))};
I_exact = [144, 2, besseli(0,2), 4];

figure;
for i=1:4
    a = A(i);
    b = B(i);
    f = F{i};
    exact_value = I_exact(i);
    error = zeros(100,1);
    
    for n=1:100
        x = linspace(a, b, n+1);
        I_trapezoid = 0;
        for j=1:n
            I_trapezoid = I_trapezoid + ((b-a)/(2*n)) * (f(x(j)) + f(x(j+1)));
        end
        error(n) = abs(I_trapezoid - exact_value);
    end
    
    % Set all zeroes to minimum value to avoid plotting empty values in place
    % of log(0).
    error(~error) = min(error(error ~= 0));
    
    loglog(linspace(1,100, 100), error, 'DisplayName', func2str(f));
    hold on; % Must hold on after plotting to keep using log scale.
end

legend('show');