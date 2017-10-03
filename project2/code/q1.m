function q1()
rho_max = 1.0;
v_max = 1.0;

function F = F(rho_L, rho_R)
    function f = f(rho)
        f = v_max * (rho - (rho^3 / rho_max^2));
    end
    
    if rho_L <= rho_max/sqrt(3) && rho_R <= rho_max/sqrt(3)
        % We include the = case here.
        F_min = f(rho_L);
        F_max = f(rho_R);
    elseif rho_L > rho_max/sqrt(3) && rho_R > rho_max/sqrt(3)
        F_min = f(rho_R);
        F_max = f(rho_L);
    elseif rho_L < rho_max/sqrt(3) && rho_max/sqrt(3) < rho_R
        F_min = min(f(rho_L), f(rho_R));
    elseif rho_L > rho_max/sqrt(3) && rho_max/sqrt(3) < rho_R
        F_max = f(rho_max/sqrt(3));
    else
        F_max = max(f(rho_L), f(rho_R));
        F_min = max(f(rho_L), f(rho_R));
    end

    % TODO: what if you put all the above cases nested under here? That way
    % you can assume assume which is larger than the other.
    if rho_L <= rho_R
        F = F_min;
    elseif rho_L > rho_R
        F = F_max;
    end
end

N = 500;
x = linspace(0,10,N+1);
dx = x(2:end) - x(1:end-1);
dt = 0.01;
n_step = 2/dt;

% initial condition
rho_0 = 0.8*rho_max*ones(N+1,1); % light traffic
rho = rho_0;

for i=1:n_step
    plot(x,rho);
    title(i*dt)
    drawnow;
    rho(1) = rho_0(1);
    
    F_imh = zeros(N+1,1);
    F_iph = zeros(N+1,1);
    for j=2:N
        F_imh(j) = F(rho(j-1), rho(j));
        F_iph(j) = F(rho(j), rho(j+1));
%         if i*dt < 2
%             F_imh(int64(N/2)) = 0;
%             F_iph(int64(N/2)) = 0;
%         end
    end
    
    if i*dt < 1
        F_imh(int64(N/2)+1) = 0;
        F_iph(int64(N/2)) = 0;
    end
    
    rho = rho - (dt./dx) .* (F_iph - F_imh);
end

end