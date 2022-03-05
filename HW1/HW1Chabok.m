clear
clc
%% physical properties

% Temperaturea based on degrees Celsius 
T_0 = 25;
T_b = 200;
% Dimensions of fin
w = 0.3;
L = 0.1;
t = 0.005;
% Heat transfer coefficients
h = 50;
k = 300;
% Fin surround and area
A_c = w * t;
P = 2 *(w + t);

%% Change of the variable

teta = T_b - T_0;
m = sqrt( h * P / (k * A_c));

%% numerical solution

% Domain steps
dx = [0.1, 0.025, 0.01, 0.001, 0.0001]; 
d2x = dx.*dx;

%Building matrix and ploting Temperature distribution
for i = 1:length(dx) 
    
   x = 0:dx(i):L; % Domain x
   n = length(x); % Number of domain segments
   % Specify the size of diameter elements of the matrix
   u = zeros(1,n); 
   d = zeros(1,n); 
   l = zeros(1,n); 
   b = zeros(1,n); 
   % Specify some matrix elements
   u(1) = 0; u(n) = 0; 
   l(1) = 1; l(n) = -1; 
   d(1) = 1; d(n) = (h/k)*dx(i)+1; 
   
          for j = 2:n-1
                l(j) = 1;
                d(j) = -m^2 * d2x(i)-2;
                u(j) = 1;
          end
     
    b(1) = T_b; % Specify B matrix element
    A = matrixillustrator(l,d,u); 
    tdma = TDMAsolver(A,b);
    s = char('b ','m ','c','k','r'); % plotting color
    plot(x,tdma,s(i,:),'linewidth',1) % plotting
            hold on
            title('Temperature distribution of numerical solution')
            axis([0,0.1,140,230])
            legendInfo{i} = ['h =',num2str(dx(i)),'m'];
            legend(legendInfo)
            xlabel('x_f')
            ylabel('T_{\circC}')
                                           
end

%% Analytical solution

x_analytical = 0:0.0001:L; % Domain x
% Analytical temperature distribution constant:
C_1 = -teta * (h * cosh(m * L) + k * m * sinh(m * L))/( k * m * cosh(m * L) + h * sinh( m * L));
% Analytical temperature distribution constant:
C_2 = teta;
% Analytical temperature distribution:
F_Teta = C_1 * sinh(m * x_analytical) + C_2 * cosh( m * x_analytical); 
% Plotting
figure (2)
plot(x_analytical,F_Teta + T_0,"b--")
title('Temperature distribution of analytical solution')
axis([0,0.1,140,230])
xlabel('x_f')
ylabel('T_{\circC}')


