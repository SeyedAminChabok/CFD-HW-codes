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

%% Successive error  in x=L with considering analytical solution and without it

x = L;

% Analytical temperature distribution constant:
C_1 = -teta * (h * cosh(m * L) + k * m * sinh(m * L))/( k * m * cosh(m * L) + h * sinh( m * L));
% Analytical temperature distribution constant:
C_2 = teta; 
% Analytical temperature distribution:
F_Teta = C_1 * sinh(m * x) + C_2 * cosh( m * x); 
% Second derivative of the analytical temperature distribution:
F2_Teta = C_1 * (m^2) * sinh(m * x) + C_2 * (m^2) * cosh( m * x);

% Domain steps
dx = [0.1, 0.025, 0.01, 0.001, 0.0001];

for i=1:length(dx)
           
     % The central finite difference of the second derivative:
     F2_Teta_numerical(i)= ((C_1*sinh(m*(x-dx(i)))+C_2*cosh(m*(x-dx(i))))-2*(C_1*sinh(m*x)+C_2*cosh(m*x))+C_1*sinh(m*(x+dx(i)))+C_2*cosh(m*(x+dx(i))))/((dx(i))^2);
     % Successive error with considering analytical solution:
     Successive_error1 = abs(F2_Teta-F2_Teta_numerical);
     % Successive error without considering analytical solution:
     Successive_error2 = - diff (F2_Teta_numerical);
     b = - diff (dx);
     
end

%% plotting

figure (1);
loglog(dx,Successive_error1,'color','r');
title('Successive error with considering analytical solution')
xlabel('h');
ylabel('Successive error');

figure (2);
loglog(b,Successive_error2,'color','b');
title('Successive error without considering analytical solution')
xlabel('h');
ylabel('Successive error');



