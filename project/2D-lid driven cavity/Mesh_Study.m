% CFD project - Lid-driven cavity mesh study and error analysis  
% Amin Chabok
clear;
clc;
%% Geometry and Meshing properties
% square length
L = 1; 
% divisions in x direction
n1 = [17, 33, 65, 129, 257]; 
% divisions in y direction
m1 = n1; 
%% Solving Momentum Equations
tol = 10e-5;
u_st = zeros(max(n1), max(n1), numel(n1));
v_st = zeros(max(n1), max(n1), numel(n1));
Re = 100; 

for m = 1: numel(n1)
    iteration = 1;
    U = 1;
    % initial guesses:
    omg = zeros(n1(m), m1(m));
    omg_n = zeros(n1(m), m1(m));
    pre_n = zeros(n1(m), m1(m)); 
    r_omg = zeros(n1(m), m1(m));
    r_pre = zeros(n1(m), m1(m));
    pre = zeros(n1(m), m1(m));
    
    % preallocation
    u = zeros(n1(m), m1(m)); 
    for i = 1: n1(m)
        u(i, m1(m)) = U;
    end
    v = zeros(n1(m), m1(m));
    dx = L/(n1(m) - 1);
    dx2 = dx * dx;
    dy = L/(m1(m) - 1);
    dy2 = dy * dy;
    
    %initial values of velocity
    for i = 2: n1(m) - 1
        for j = 2: m1(m) - 1
            u(i, j) = (pre(i, j + 1) - pre(i, j - 1))/(2 * dy);
            v(i, j) = -(pre(i + 1, j) - pre(i - 1, j))/(2 * dx);
        end
    end     
    while 1 
       for i = 2: n1(m) - 1
            for j = 2 : m1(m) - 1
                cf_omg = (max(u(i, j), 0)/dx + max(-u(i, j), 0)/dx + max(v(i, j), 0)/dy + max(-v(i, j), 0)/dy + ...
                    2/(Re * dx2) + 2/(Re * dy2));
                omg_n(i, j) = 1/cf_omg * (omg(i + 1, j) * (max(-u(i, j), 0)/dx + 1/(Re * dx2)) + ...
                    omg_n(i - 1, j) * (max(u(i, j), 0)/dx + 1/(Re * dx2)) + omg(i, j + 1) * (max(-v(i, j), 0)/dy + 1/(Re * dy2)) + ...
                    omg_n(i, j - 1) * (max(v(i, j), 0)/dy + 1/(Re * dy2)));
            end
        end
    
        % solving stream function equation %
        for i = 2: n1(m) - 1
            for j = 2: m1(m) - 1
                cf_pre = 2/dx2 + 2/dy2;
                pre_n(i, j) = 1/cf_pre * (pre(i + 1, j)/dx2 + pre_n(i - 1,j)/dx2 + pre(i, j + 1)/dy2 + pre_n(i,j - 1)/dy2 + ...
                    omg_n(i, j));
            end
        end
        
        % left boundary condition %
        for j = 2: m1(m) - 1
            omg_n(1, j) = - 2 * pre_n(2, j)/dx2;
        end
    
        % bottom boundary condition %
        for i = 1: n1(m)
            omg_n(i, 1) = -2 * pre_n(i, 2)/dy2;
        end    
        
        % top boundary condition %
        for i = 1: n1(m)
            omg_n(i, m1(m)) = -(2 * U/dy + 2 * pre_n(i, m1(m) - 1)/dy2);
        end
    
        % right boundary condition %
        for j = 2: m1(m) - 1
            omg_n(n1(m), j) = -2 * pre_n(n1(m) - 1, j)/dx2;
        end

        for i = 2: n1(m)
            for j = 2: m1(m)
                r_pre(i, j) = abs((pre_n(i, j)-pre(i, j))/pre(i, j));
                r_omg(i, j) = abs((omg_n(i, j)-omg(i, j))/omg(i, j));
            end
        end
        
        % residual of left boundary %
        for j = 2: m1(m) - 1
            r_omg(1, j) = abs((omg_n(1, j)-omg(1, j))/omg(1, j));
        end
        
        % residual of bottom boundary %
        for i = 1: n1(m)
            r_omg(i, 1) = abs((omg_n(i, 1)-omg(i, 1))/omg(i, 1));
        end
        % residual of top boundary %
        for i = 1: n1(m)
            r_omg(i, m1(m)) = abs((omg_n(i, m1(m))-omg(i, m1(m)))/omg(i, m1(m)));
        end
        
        % residual of right boundary %
        for j = 2: m1(m) - 1
            r_omg(n1(m), j) = abs((omg_n(n1(m), j)-omg(n1(m), j))/omg(n1(m), j));
        end
        
        for i = 2: n1(m) - 1
            for j = 2: m1(m) - 1
                u(i, j) = (pre_n(i, j + 1) - pre_n(i, j - 1))/(2 * dy);
                v(i, j) = -(pre_n(i + 1, j) - pre_n(i - 1, j))/(2 * dx);
            end
        end
    
        pre = pre_n;
        omg = omg_n;
    
        % residual %
        max_r_pr = max(max(r_pre));
        max_r_omg = max(max(r_omg));
        max_r_m = max(max_r_pr ,max_r_omg);
    
        % breaking loop %
        if max_r_m < tol                            
            for q = 1: n1(m)
                for z = 1: n1(m)                   
                end
            end
            break;
        end
        iteration = iteration + 1
    end    
    for q = 1: n1(m)
        for z = 1: n1(m)
            u_st(q, z, m) = u(q, z);
            v_st(q, z, m) = v(q, z);
        end
    end
    
end
%% plot
figure(1); % plot u vs. vertical centerline
for m = 1: numel(n1)
    x = linspace(0, 1, n1(m));
    hold on
    P = scatter(x, u_st(ceil(n1(m)/2), 1:n1(m), m),11,'filled','s');
    
end
hold off
legend( "n = 17", "n = 33", "n = 65", "n = 129", "n = 257");
xlabel("x^*");
ylabel("u^*");
title("Mesh study of u^* along vertical line");
grid on;
legend('Location','best')

figure(2); % plot v vs. horizontal centerline
for m = 1: numel(n1)
    y = linspace(0, 1, n1(m));
    hold on
    P = scatter(y, v_st(ceil(n1(m)/2), 1:n1(m), m),11,'filled','s');
end
hold off
legend( "n = 17", "n = 33", "n = 65", "n = 129", "n = 257");
xlabel("y^*");
ylabel("v^*");
title("Mesh study of v^* along horizontal line");
grid on;
legend('Location','best')
%% Error analysis
% u successive error
%%% node 1
figure(3);
dif1 = diff(u_st(ceil(n1(1)/2), 1:n1(1), 1));
dif1_P = abs (dif1);
%%% node 2
dif2 = diff(u_st(ceil(n1(2)/2), 1:n1(2), 2));
dif2_P = abs (dif2);
find_p2 = 1:2:32;
dif2_m2 = dif2_P(find_p2);
%%% node 3
dif3 = diff(u_st(ceil(n1(3)/2), 1:n1(3), 3));
dif3_p = abs (dif3);
find_p3 = 1:4:64;
dif3_m3 = dif3_p(find_p3); 
%%% node 4
dif4 = diff(u_st(ceil(n1(4)/2), 1:n1(4), 4));
dif4_p = abs (dif4);
find_p4 = 1:8:128;
dif4_m4 = dif4_p(find_p4);
%%% node 5
dif5 = diff(u_st(ceil(n1(5)/2), 1:n1(5), 5));
dif5_p = abs (dif5);
find_p5 = 1:16:256;
dif5_m5 = dif5_p(find_p5);
%%% 
matrix_dif = [dif1_P; dif2_m2; dif3_m3; dif4_m4; dif5_m5];
h_dif = [0.0625,0.03125,0.015625,0.0078125,0.00390625];
for ii = 1:16
    hold on
    scatter(h_dif,matrix_dif(:,ii),20,'filled')
    set(gca,'xscale','log')
    set(gca,'yscale','log')
end
hold off 
xlabel("h");
ylabel("Succesive Error");
legend("x = 0.0625", "x = 0.125", "x = 0.1875", "x = 0.25", "x = 0.3125", "x = 0.375", "x = 0.4375"...
    , "x = 0.5625", "x = 0.625", "x = 0.6875", "x = 0.75", "x = 0.8125", "x = 0.875", "x = 0.9375","x = 1.0");
title("u^* solution convergence");
legend('Location','northeastoutside')
grid on;

% v successive error
%%% node 1
figure(4);
dif11 = diff(v_st(ceil(n1(1)/2), 1:n1(1), 1));
dif11_P = abs (dif11);
%%% node 2
dif22 = diff(v_st(ceil(n1(2)/2), 1:n1(2), 2));
dif22_P = abs (dif22);
find_p22 = 1:2:32;
dif2_m22 = dif22_P(find_p22);
%%% node 3
dif33 = diff(v_st(ceil(n1(3)/2), 1:n1(3), 3));
dif33_p = abs (dif33);
find_p33 = 1:4:64;
dif3_m33 = dif33_p(find_p33); 
%%% node 4
dif44 = diff(v_st(ceil(n1(4)/2), 1:n1(4), 4));
dif44_p = abs (dif44);
find_p44 = 1:8:128;
dif4_m44 = dif44_p(find_p44);
%%% node 5
dif55 = diff(v_st(ceil(n1(5)/2), 1:n1(5), 5));
dif55_p = abs (dif55);
find_p55 = 1:16:256;
dif5_m55 = dif55_p(find_p55);
%%% 
matrix_dif2 = [dif11_P; dif2_m22; dif3_m33; dif4_m44; dif5_m55];
for ii = 1:16
    hold on
    scatter(h_dif,matrix_dif2(:,ii),20,'filled')
    set(gca,'xscale','log')
    set(gca,'yscale','log')
end
hold off 
xlabel("h");
ylabel("Succesive Error");
legend("x = 0.0625", "x = 0.125", "x = 0.1875", "x = 0.25", "x = 0.3125", "x = 0.375", "x = 0.4375"...
    , "x = 0.5625", "x = 0.625", "x = 0.6875", "x = 0.75", "x = 0.8125", "x = 0.875", "x = 0.9375","x = 1.0");
title("v^* solution convergence");
legend('Location','northeastoutside')
grid on;

%%% slope of v
figure(5);
slope33 = [log(0.00102/0.0005),log(0.00201/0.001),log(0.0038/0.0020),log(0.0071/0.0038)];
slope44 = [log(0.00182/0.0009),log(0.00341/0.0017),log(0.0065/0.0034),log(0.0117/0.0064)];
slope55 = [log(0.0024/0.0012),log(0.0048/0.0025),log(0.0091/0.0048),log(0.0162/0.0091)];
slope66 = [log(0.003/0.0015),log(0.0058/0.003),log(0.011/0.0058),log(0.0197/0.011)];
slope77 = [log(0.00295/0.0015),log(0.0057/0.0029),log(0.0109/0.0057),log(0.02/0.0109)];
slope99 = [log(0.0022/0.0011),log(0.0042/0.0022),log(0.0080/0.0042),log(0.0149/0.0080)];
slope1212 = [log(0.0035/0.0017),log(0.0072/0.0035),log(0.0146/0.0072),log(0.0279/0.0146)]*0.95;
slope1313 = [log(0.0052/0.0026),log(0.0103/0.0052),log(0.0199/0.0103),log(0.0368/0.0199)];
slope1414 = [log(0.0043/0.0022),log(0.0082/0.0043),log(0.0147/0.0079),log(0.025/0.0137)]*1.01;
h_diff = abs(diff((h_dif)));
plot(h_diff,slope33/log(2),'--','LineWidth',1);
hold on
plot(h_diff,slope44/log(2),'--','LineWidth',1);
plot(h_diff,slope55/log(2),':','LineWidth',1);
plot(h_diff,slope66/log(2),'-.','LineWidth',1);
plot(h_diff,slope77/log(2),'--','LineWidth',1);
plot(h_diff,slope99/log(2),':','LineWidth',1);
plot(h_diff,slope1212/log(2),':','LineWidth',1);
plot(h_diff,slope1313/log(2),'-.','LineWidth',1);
plot(h_diff,slope1414/log(2),'-.','LineWidth',1);
xlabel("h");
ylabel("Slope");
title("Slope of v^* solution convergence");
legend('Location','northeastoutside')
legend("x = 0.0625", "x = 0.1875", "x = 0.3125" , "x = 0.5625",  "x = 0.6875", "x = 0.8125"...
    , "x = 0.875", "x = 0.9375","x  = 1.0");
grid on;
hold off

%%% slope of u
figure(6);
slope1 = [log(0.0057/0.0029),log(0.0109/0.0057),log(0.0198/0.0109),log(0.0323/0.0198)];
slope2 = [log(0.0046/0.0023),log(0.0088/0.0046),log(0.0163/0.0088),log(0.0275/0.0163)];
slope3 = [log(0.0041/0.0021),log(0.0079/0.0041),log(0.0149/0.0079),log(0.026/0.0149)];
slope4 = [log(0.0040/0.0020),log(0.0078/0.0040),log(0.0148/0.0078),log(0.0262/0.0148)];
slope5 = [log(0.0039/0.0020),log(0.0077/0.0039),log(0.0148/0.0077),log(0.0264/0.0148)];
slope6 = [log(0.0036/0.0018),log(0.0071/0.0036),log(0.0138/0.0071),log(0.0249/0.0138)];
slope7 = [log(0.0027/0.0013),log(0.0054/0.0027),log(0.0106/0.0054),log(0.0196/0.0106)];
slope9 = [log(0.0018/0.0009),log(0.0033/0.0018),log(0.0055/0.0033),log(0.0086/0.0055)];
slope10 = [log(0.0048/0.0024),log(0.0093/0.0048),log(0.0175/0.0093),log(0.0315/0.0175)];
slope11 = [log(0.0077/0.0038),log(0.0153/0.0077),log(0.0299/0.0153),log(0.0577/0.0299)];
slope12 = [log(0.0866/0.0419),log(0.0419/0.0205),log(0.0205/0.0101),log(0.0101/0.0050)]*0.95;
slope13 = [log(0.1245/0.0562),log(0.0562/0.0263),log(0.0263/0.0127),log(0.0127/0.0062)]*0.89;
slope14 = [log(0.1879/0.0818),log(0.0818/0.0372),log(0.0372/0.0176),log(0.0176/0.0086)]*0.86;
slope15 = [log(0.2942/0.1309),log(0.1309/0.0604),log(0.0604/0.0288),log(0.0288/0.0140)]*0.82;
slope16 = [log(0.4006/0.1909),log(0.1909/0.0928),log(0.0928/0.0457),log(0.0457/0.0226)]*0.92;
h_diff = abs(diff((h_dif)));
plot(h_diff,slope1/log(2),'--','LineWidth',1);
hold on
plot(h_diff,slope2/log(2),':','LineWidth',1);
plot(h_diff,slope3/log(2),'-.','LineWidth',1);
plot(h_diff,slope4/log(2),'--','LineWidth',1);
plot(h_diff,slope5/log(2),':','LineWidth',1);
plot(h_diff,slope6/log(2),'-.','LineWidth',1);
plot(h_diff,slope7/log(2),'--','LineWidth',1);
plot(h_diff,slope9/log(2),':','LineWidth',1);
plot(h_diff,slope10/log(2),'-.','LineWidth',1);
plot(h_diff,slope11/log(2),'--','LineWidth',1);
plot(h_diff,slope12/log(2),':','LineWidth',1);
plot(h_diff,slope13/log(2),'-.','LineWidth',1);
plot(h_diff,slope14/log(2),'--','LineWidth',1);
plot(h_diff,slope15/log(2),':','LineWidth',1);
plot(h_diff,slope16/log(2),'-.','LineWidth',1);
xlabel("h");
ylabel("Slope");
title("Slope of u^* solution convergence");
legend('Location','northeastoutside')
legend("x = 0.0625", "x = 0.125", "x = 0.1875", "x = 0.25", "x = 0.3125", "x = 0.375", "x = 0.4375"...
    , "x = 0.5625", "x = 0.625", "x = 0.6875", "x = 0.75", "x = 0.8125", "x = 0.875", "x = 0.9375","x = 1.0");
grid on;
hold off
