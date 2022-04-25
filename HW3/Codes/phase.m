%% phase
L = 1.0 ;
dx = L/50.0 ;
x = 0 : dx :  L ;
nx = length(x) ;
cfl = [0.0025,0.0050,0.00625,0.0125,1] ;

for icfl = 1 : length(cfl)
  for ix = 1 : nx
            k(ix) = 2.0*pi*ix/(2.0*L) ;
            beta(ix) = k(ix)*dx ;
            amplitude(ix) = cos(beta(ix))^2.0 + cfl(icfl)*sin(beta(ix))^2.0 ;
         %  phi(ix) = atan(-cfl(icfl)*tan(beta(ix))) ;
      end
    
    figure(2)
    hold on
    grid on
    plot (x,amplitude,'-s','linewidth',1.2)
    xlabel 'relative phase'
    ylabel amplitude
end
legend('CFL = 0.0025','CFL = 0.0050','CFL = 0.00625','CFL = 0.0125','CFL = 1.0') 
