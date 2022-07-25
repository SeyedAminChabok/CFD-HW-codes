% CFD project - Lid-driven cavity validation  
% Amin Chabok
clear;
clc;
%% Geometry and Meshing properties
L = 1; 
% divisions in x direction
N = 129;
% divisions in y direction
M = 129;
x = linspace(0, L, N);
dx = L/(N - 1);
dx2 = dx * dx;
y = linspace(0, L, M);
dy = L/(M - 1);
dy2 = dy * dy;
%% Solving Momentums Equations
tol = 10e-6;
Re = 100;
U = 1; 
iteration_m = 1;
% initial guesses:
omg = zeros(N, M); 
omg_n = zeros(N, M);
pr = zeros(N, M); 
pr_n = zeros(N, M);
r_pr = zeros(N, M);
r_omg = zeros(N, M); 
u = zeros(N, M); 
for i = 1: N
    u(i, M) = U;
end
v = zeros(N, M);  

% preallocation
for i = 2: N - 1
    for j = 2: M - 1 
        u(i, j) = (pr(i, j + 1) - pr(i, j - 1))/(2 * dy);
        v(i, j) = -(pr(i + 1, j) - pr(i - 1, j))/(2 * dx);
    end
end

while 1    
    for i = 2: N - 1
        for j = 2 :M - 1
           cf_omg = (max(u(i, j), 0)/dx + max(-u(i, j), 0)/dx + max(v(i, j), 0)/dy + max(-v(i, j), 0)/dy + ... 
               2/(Re * dx2) + 2/(Re * dy2));
           omg_n(i, j) = 1/cf_omg * (omg(i + 1, j) * (max(-u(i, j), 0)/dx + 1/(Re * dx2)) + ...
               omg_n(i - 1, j) * (max(u(i, j), 0)/dx + 1/(Re * dx2)) + omg(i, j + 1) * (max(-v(i, j), 0)/dy + 1/(Re * dy2)) + ...
               omg_n(i, j - 1) * (max(v(i, j), 0)/dy + 1/(Re * dy2)));
        end
    end
    
    for i = 2: N - 1
        for j = 2: M - 1
            cf_pr = 2/dx2 + 2/dy2;
            pr_n(i, j) = 1/cf_pr * (pr(i + 1, j)/dx2 + pr_n(i - 1,j)/dx2 + pr(i, j + 1)/dy2 + pr_n(i,j - 1)/dy2 + ...
                omg_n(i, j));
        end
    end
    
     % left boundary: i = 1, j = 2: M - 1
    for j = 2: M - 1
        omg_n(1, j) = - 2 * pr_n(2, j)/dx2;
    end
    
    % bottom boundary: j = 1, i = 1: N
    for i = 1: N
        omg_n(i, 1) = -2 * pr_n(i, 2)/dy2;
    end   
    % top boundary: j = M, i = 1: N
    for i = 1: N
        omg_n(i, M) = -(2 * U/dy + 2 * pr_n(i, M - 1)/dy2);
    end
    
    % right boundary: i = N, j = 2: M - 1
    for j = 2: M - 1
        omg_n(N, j) = -2 * pr_n(N - 1, j)/dx2;
    end

    for i = 2: N - 1
        for j = 2: M - 1
            r_pr(i, j) = abs((pr_n(i, j)-pr(i, j))/pr(i, j));
            r_omg(i, j) = abs((omg_n(i, j)-omg(i, j))/omg(i, j));
        end
    end
    
    % residual of left boundary
    for j = 2: M - 1
        r_omg(1, j) = abs((omg_n(1, j)-omg(1, j))/omg(1, j));
    end
    
    % residual of bottom boundary
    for i = 1: N
        r_omg(i, 1) = abs((omg_n(i, 1)-omg(i, 1))/omg(i, 1));
    end
    % residual of top boundary
    for i = 1: N
        r_omg(i, M) = abs((omg_n(i, M)-omg(i, M))/omg(i, M));
    end
    
    % residual of right boundary 
    for j = 2: M - 1
        r_omg(N, j) = abs((omg_n(N, j)-omg(N, j))/omg(N, j));
    end
    
    for i = 2: N - 1
        for j = 2: M - 1
            u(i, j) = (pr_n(i, j + 1) - pr_n(i, j - 1))/(2 * dy);
            v(i, j) = -(pr_n(i + 1, j) - pr_n(i - 1, j))/(2 * dx);
        end
    end
    pr = pr_n;
    omg = omg_n;
    max_r_pr = max(max(r_pr));
    max_r_omg = max(max(r_omg));
    max_r_m = max(max_r_pr ,max_r_omg);

    if max_r_m < tol
        break;
    end
    iteration_m = iteration_m + 1
end

%% Ghia paper data
u_g = [0, -0.03717, -.04192, -.04775, -.06434, -.10150, -.15662, -.21090, -.20581, -.13641, .00332, .23151, .68717, .73722, ... 
    .7887, .84123, 1];
v_g = [0, .09233, .10091, .10890, .12317, .16077, .17507, .17527, .05454, -.24533, -.22445, -.16914, -.10313, -.08864, -.07391,... 
    -.05906, 0];
u_s = [u(65, 1), u(65, 8), u(65, 9), u(65, 10), u(65, 14), u(65, 23), u(65, 37), u(65, 59), u(65, 65), u(65, 80),...
    u(65, 95), u(65, 110), u(65, 123), u(65, 124), u(65, 125), u(65, 126), u(65, 129)];
v_s = [v(1, 65), v(9, 65), v(10, 65), v(11, 65), v(13, 65), v(21, 65), v(30, 65), v(31, 65), v(65, 65), v(104, 65),...
    v(111, 65), v(117, 65), v(122, 65), v(123, 65), v(124, 65), v(125, 65), v(129, 65)];
x_s = [x(1), x(9), x(10), x(11), x(13), x(21), x(30), x(31), x(65), x(104), x(111), x(117), x(122), x(123), x(124), x(125), x(129)];
y_s = [y(1), y(8), y(9), y(10), y(14), y(23), y(37), y(59), y(65), y(80), y(95), y(110), y(123), y(124), y(125), y(126), y(129)];

%% plots
figure(1)
plot(y_s, u_g, 'ro', y_s, u_s, 'b');
grid on;
xlabel("y");
ylabel("u");
title("Validation for u-velocity along Vertical Line through Geometric Center of Cavity - Re = 100");
legend("Ghia", "Current study")
legend('Location','best')

figure(2)
plot(x_s, v_g, 'ro', x_s, v_s, 'b');
grid on;
xlabel("x");
ylabel("v");
title("Validation for v-Velocity along Horizontal Line through Geometric Center of Cavity - Re = 100");
legend("Ghia", "Current study")
legend('Location','best')


