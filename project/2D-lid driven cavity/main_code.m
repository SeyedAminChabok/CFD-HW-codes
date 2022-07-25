% CFD project - Lid-driven cavity main code
% Amin Chabok
clear;
clc;
%% Geometry and Meshing properties
L = 1;
Re = 10;  % Re = 1,10,100,500 
N = 129;
M = 129;
U=1;
dx=1/(N-1);
dy=1/(M-1);
factor=2*(1/dx^2+1/dy^2);
re=1;
reb=re;
%% Solving Momentum Equations
p = zeros(N,M);
o = zeros(N,M);
ss = zeros(N,M);       
w = zeros(N,M);         
uac = zeros(N,M);
vac = zeros(N,M);
u = zeros(N,M);
v = zeros(N,M);
u(:,M) = 1;
pr = zeros(N,M);        
prr = zeros(N,M);
tol = 1;
iteration = 0;
    
    while tol>10^-5
        tol=0;
        iteration = iteration + 1
        
        for i=2:N-1
            for j=2:M-1                
                pold(i,j)=p(i,j);
                fj=1/factor*(o(i,j)+(p(i+1,j)+p(i-1,j))/dx^2+(p(i,j+1) ...
                    +p(i,j-1))/dy^2);
                p(i,j)=p(i,j)+re*(fj-p(i,j));
                
            end
        end
        
        for j=1:M
            ads=-2*p(2,j)/dx^2;     
            o(1,j)=o(1,j)+reb*(ads-o(1,j));
            fcd=-2*p(N-1,j)/dx^2;  
            o(N,j)=o(N,j)+reb*(fcd-o(N,j));
        end
        for i=1:N
            mnu=-2*p(i,2)/dy^2;                 
            o(i,1)=o(i,1)+reb*(mnu-o(i,1));
            mnn=-(2*p(i,M-1)+2*U*dy)/dy^2;     
            o(i,M)=o(i,M)+reb*(mnn-o(i,M));
        end
        for i=2:N-1
            for j=2:M-1
                od(i,j)=o(i,j);
                fh=1/factor*((o(i+1,j)+o(i-1,j))/dx^2+(o(i,j+1)+o(i,j-1))/dy^2 ...
                -Re*(p(i,j+1)-p(i,j-1))*(o(i+1,j)-o(i-1,j))/(4*dx*dy) ...
                +Re*(p(i+1,j)-p(i-1,j))*(o(i,j+1)-o(i,j-1))/(4*dx*dy));
                
                o(i,j)=o(i,j)+re*(fh-o(i,j));
                tol=tol+abs(o(i,j)-od(i,j));
                
            end
        end
        
   
        for i=2:N-1
            for j=2:M-1
                u(i,j)=(p(i,j+1)-p(i,j-1))/(2*dy);
                v(i,j)=-(p(i+1,j)-p(i-1,j))/(2*dx);
            end
        end
     
     
        for i=(2:N-1)
            for j=(2:M-1)
                prr(i,j)=0.25*((pr(i+1,j)+pr(i-1,j)+pr(i,j+1)+pr(i,j-1)...
                    -(0.5*((u(i+1,j)-u(i-1,j))*((v(i,j+1)-v(i,j-1))...
                    -((u(i,j+1)-u(i,j-1))*((v(i+1,j)-v(i-1,j)))))))));                
            end
        end
        
        pr=prr;
        
        for i=1:N
            for j=1:M
                k=i;    l=j;
                ss(l,k)=p(i,j);
                w(l,k)=o(i,j);
                uac(l,k)=u(i,j);
                vac(l,k)=v(i,j);
                prc(l,k)=pr(i,j);
            end
        end
    end
%% plots

figure(1)    
Z=ss(1:N,1:M);  
L=linspace(0,1,N);
L=linspace(0,1,M);
[c,h]= contour(L,L,Z);
xlabel("x^*");
ylabel("y^*");
title('stream function contours - Re=10')
colorbar
    
figure(2)  
[m,n]=contour(L,L,prc,100);
xlabel("x^*");
ylabel("y^*");
title('Pressure contours - Re=10')
colorbar

figure(3)  
[m,n]=contour(L,L,uac,100);
xlabel("x^*");
ylabel("y^*");;
title('u distribution - Re=10')
colorbar
    

figure(4)  
[m,n]=contour(L,L,vac,100);
xlabel("x^*");
ylabel("y^*");
title('v distribution - Re=10')
colorbar
