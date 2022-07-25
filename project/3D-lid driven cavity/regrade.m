function y=regrade(x,a,e)
%REGRADE  grade nodal points in array towards edge or center
%
% Regrades array of points
%
% Usage: regrade(x,a,e)
% x is array of nodal point coordinates in increasing order.
% a is parameter which controls grading.
% e selects side or sides for refinement.
%
% if e=0: refine both sides, 1: refine upper, 2: refine lower.
%
% if a=1 then return xarray unaltered.
% if a<1 then grade towards the edge(s)
% if a>1 then grade away from edge.

ae=abs(a);
n=length(x);
y=x;
if ae==1 | n<3 | (e~=0 & e~=1 & e~=2) return;
end;
if e==0 
  Xmx=max(x);
  Xmn=min(x);
  Xc=(Xmx+Xmn)/2;
  Xl=(Xmx-Xmn)/2;
  for k=2:(n-1)
    xk=x(k)-Xc;
    y(k)=Xc+Xl*sign(xk)*(abs(xk/Xl))^ae;
 end
elseif (e==1 & x(1)<x(n)) | (e==2 & x(1)>x(n))
   for k=2:n-1
      xk=x(k)-x(1);
      y(k)=x(1)+(x(n)-x(1))*(abs(xk/(x(n)-x(1))))^ae;
   end
else   % (e==2 & x(1)<x(n)) | (e==1 & x(1)>x(n))
   for k=2:n-1
      xk=x(k)-x(n);
      y(k)=x(n)+(x(1)-x(n))*(abs(xk/(x(1)-x(n))))^ae;
   end
end
return;