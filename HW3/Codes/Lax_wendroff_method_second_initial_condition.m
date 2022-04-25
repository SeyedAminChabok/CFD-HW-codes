clc, clear, close all ;
%% second initial condition 

dts=[0.00050,0.00025,0.00020,0.00010];
for i = 1 : length(dts)
   t = 0.0 ; 
  dt = dts(i) ;
  nt = 1000 ;
    
  L = 1.0 ;
  dx = L/50.0 ;
  dxx = 2.0*dx ;
  ddx = dx*dx ;
  a = 0.5 ;
  
  x = 0.0 : dx : L ;
  nx = length(x) ;
  
  u = zeros(nx,1) ;
  u0 = zeros(nx,1) ;
  
  CFL = [a*dts(1)/dx, a*dts(2)/dx, a*dts(3)/dx, a*dts(4)/dx] ;
  
  %% initial condition
  for ix = 1 : nx
        u0(ix) = sin(4.0*pi*x(ix)) ;
  end
  
  % PBC 
  for ix = 1 : nx
        ixp(ix) = ix+1 ;
        ixm(ix) = ix-1 ;
  end
    ixp(nx) = 2 ;
    ixm(1) = nx-1 ;
  
  for it = 1 : nt
      for ix = 1 : nx
        u(ix) = u0(ix)-a*dt*((u0(ixp(ix))-u0(ixm(ix)))/dxx) + 0.5*a^2.0*dt^2.0*...
                ((u0(ixp(ix))-2.0*u0(ix)+u0(ixm(ix))))/ddx ;
      end
      u0 = u ;
      
    for ix = 1 : nx
        u_exact(ix) = sin(4.0*pi*(x(ix)-a*t)) ;
    end
      t = t + dt ;
  end
  
    hold on
    grid on
    plot (x,u,'-s',x,u_exact, '-.s','linewidth',1.5)
    axis ([0.0 1.0 -1.0 1.0])
    xlabel x
    ylabel u
    title 'Lax-Wendroff Method (second initial condition)'

end
    legend('CFL = 0.0125','Exact-CFL','CFL = 0.00625',...
        'Exact-CFL = 0.00625','CFL = 0.005','Exact-CFL = 0.005','CFL = 0.0025','Exact-CFL = 0.0025')
  for ix = 1 : nx
      error(ix) = abs(u(ix)-u_exact(ix)) ;
     fprintf('\n x: %f  Numerical u: %f  Analytical u: %f Error: %f \n', x(ix), u(ix), u_exact(ix), error(ix))
  end

