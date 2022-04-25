%% phase
L = 1.0 ;
dx = L/50.0 ;
x = 0.0 : dx :  L ;
nx = length(x) ;
cfl = [0.25,0.50,0.75,1] ;

for icfl = 1 : length(cfl)
  for ix = 1 : nx
            k(ix) = 2.0*pi*ix/(2.0*L) ;
            beta(ix) = k(ix)*dx ;
            amplitude(ix) = ((cos(beta(ix)))^2.0 + cfl(icfl)*(sin(beta(ix)))^2.0)^0.5 ;
            phi(ix) = atan(-cfl(icfl)*tan(beta(ix))) ;
            a(ix) = cfl(icfl)*sin(beta(ix)) ;
            b(ix) = cos(beta(ix));
      end
    
    figure(2)
    hold on
    grid on
    plot (b,a,'-s','linewidth',1.2)
    xlabel 'relative phase'
    ylabel amplitude
end
legend('CFL = 0.25','CFL = 0.50','CFL = 0.75','CFL = 1.0') 
