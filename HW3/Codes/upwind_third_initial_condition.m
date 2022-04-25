clc, clear, close all ;

%% Third initial conditions 

dts=[0.00050,0.00025,0.00020,0.00010];
for i = 1 : length(dts)
    t = 0.0 ;
  dt = dts(i) ;
  nt = 1000 ;
    
  L = 1.0 ;
  dx = L/50.0 ;
  a = 0.5 ;
  
  x = 0.0 : dx : L ;
  nx = length(x) ;
  
  u = zeros(nx,1) ;
  u_exact = zeros(nx,1) ;
  error = zeros(nx,1) ;
  
   CFL = [a*dts(1)/dx, a*dts(2)/dx, a*dts(3)/dx, a*dts(4)/dx] ;  % CFL number 
  
  %% initial condition
  for ix = 1 : nx
      if (x(ix)>=0.0 && x(ix)<0.2)
          u(ix) = 0.0 ;
      elseif (x(ix)>=0.2 && x(ix)<0.3)
          u(ix) = 1.0 ;
      else
          u(ix) = 0.0 ;      
      end
  end
  
  RHS = 0.0*u ;
  for it = 1 : nt
      for ix = 2 : nx-1
          dudx = (u(ix)-u(ix-1))/(dx) ; % backward finite difference 
          RHS(ix) = -a*dudx ;
      end
      u = u + dt*RHS ; 
      
  for ix = 1 : nx
      if (x(ix)>=0.0 && x(ix)<0.2)
          u_exact(ix) = 0.0 ;
      elseif (x(ix)>=0.2 && x(ix)<0.3)
          u_exact(ix) = 1.0 ;
      else
          u_exact(ix) = 0.0 ;      
      end
  end
  
  end
  
    hold on
    grid on
    plot (x,u,'-s','linewidth',1.5)
    xlabel x
    ylabel u
    title 'Upwind Method (third initial condition)'
  
end
plot(x,u_exact,'-.K','linewidth',1.5)
legend('CFL = 0.0125','CFL = 0.00625','CFL = 0.005','CFL = 0.0025', 'Exact') 

   for ix = 1 : nx
       error(ix) = abs(u(ix)-u_exact(ix)) ;
     fprintf('\n x: %f  Numerical u: %f  Analytical u: %f Error: %f \n', x(ix), u(ix), u_exact(ix), error(ix))
   end

%% End of third initial condition