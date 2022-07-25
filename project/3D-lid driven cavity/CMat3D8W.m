function [Cm,RowNdx,ColNdx]=CMat3D8W(Xe, Elcon, nn2nft, Vdof)
% CMat3D8W - Returns the element convection matrix for the  
%   3D linear-complete, normal-conforming, divergence-free, Hermite basis 
%   functions on 8-node rectangular hexahedral elements with 6 DOF per node 
%   using Gauss quadrature on the 2x2x2 reference cube. 
% The 8 columns of the array V3dof must contain the six degree-of-freedom 
%   vectors in the nodal order (A,B,C,u,v,w). 
% The assumed nodal numbering starts with 1 at the lower left corner (-1,-1,-1) 
%   of the element. 
%
% Usage:
%   [Cm,Rndx,Cndx] = CMat3D8W(Xe, Elcon, nn2nft,V3dof)
%   [Cm,Rndx,Cndx,Rcm,RcNdx] = CMat3D8W(Xe, Elcon, nn2nft,V3dof)
%   Xe(1,:) -  x-coordinates of 8 corner nodes of element.  
%   Xe(2,:) -  y-coordinates of 8 corner nodes of element.  
%   Xe(3,:) -  z-coordinates of 8 corner nodes of element.  
%   Elcon(8) - connectivity matrix for this element, list of nodes. 
%   nn2nft(1,n) - global freedom number for node n.
%   nn2nft(2,n) - global freedom type for node n.
%   Vdof(6,8)  - VP & velocity Dofs at 8 nodes.
%
% Calls:
%   V8cW(nc,x,y,z), V8xyzcW(nc,x,y,z)
%
% Jonas Holdeman, July 2011
%

% Constants and fixed data
nc=[-1,-1,-1; 1,-1,-1; -1,1,-1; 1,1,-1; -1,-1,1; 1,-1,1; -1,1,1; 1,1,1]; % defines corner nodal order
ndfn=6;        % number of degrees of freedom per node. 
nne=8;         % number of nodes per element
ndfe=ndfn*nne;  % number of degrees of freedom per element

% Define 5-point quadrature data once, on first call. 
% Gaussian weights and absissas to integrate 9th degree polynomials exactly. 
global GQC5;
if (isempty(GQC5))   % Has 5-point quadrature data been defined? If not, define arguments & weights. 
   Aq=[-.906179845938664,-.538469310105683, .0,               .538469310105683, .906179845938664];
   Hq=[ .236926885056189, .478628670499366, .568888888888889, .478628670499366, .236926885056189];
   GQC5.xa=zeros(125,1); GQC5.ya=zeros(125,1); GQC5.za=zeros(125,1); GQC5.wt=zeros(125,1); 
   nr=0; 
   for nz=1:5; for ny=1:5; for nx=1:5
       nr=nr+1; GQC5.xa(nr)=Aq(nx); GQC5.ya(nr)=Aq(ny); GQC5.za(nr)=Aq(nz); 
       GQC5.wt(nr)=Hq(nx)*Hq(ny)*Hq(nz);
   end; end; end
   GQC5.size=nr; 
end

xa=GQC5.xa; ya=GQC5.ya; za=GQC5.za; W=GQC5.wt; Nq=GQC5.size;

% ---------------------------------------------------
global Z3_V8c; global Z3_V8xc; global Z3_V8yc; global Z3_V8zc;
if (isempty(Z3_V8c) | isempty(Z3_V8xc) | size(Z3_V8xc,2)~=Nq)
  Z3_V8c=cell(nne,Nq); Z3_V8xc=cell(nne,Nq); 
  Z3_V8yc=cell(nne,Nq); Z3_V8zc=cell(nne,Nq);
  for k=1:Nq
    for m=1:nne
      ncm=nc(m,:);
      Z3_V8c{m,k}=V8cW(ncm,xa(k),ya(k),za(k));
      [Z3_V8xc{m,k},Z3_V8yc{m,k},Z3_V8zc{m,k}]=V8xyzcW(ncm,xa(k),ya(k),za(k));
    end
  end
end   % if (isempty(*))
% ----------------- End fixed data ------------------
 
Ti=cell(nne);
for m=1:nne
  Jt=Xe*GTL(nc(:,:),nc(m,1),nc(m,2),nc(m,3));
  Det=det(Jt);
  JtiD=inv(Jt)*Det;
  J=Jt';
  Ti{m}=blkdiag(J,JtiD);
end  % loop m

Cm=zeros(ndfe,ndfe); S=zeros(3,ndfe);   % Preallocate arrays

% Begin loop over Gauss-Legendre quadrature points. 
for k=1:Nq  
   
  Jt=Xe*GTL(nc(:,:),xa(k),ya(k),za(k));    % transpose of Jacobian at (xa,ya,za)
  Det=det(Jt);
  JtbD=Jt/Det;
  Jti=inv(Jt); 
  JtiD=Jti*Det; 
  Ji=Jti';
%
% Compute mapped element Si and the fluid velocity at the quadrature point (xa,ya,za).
  Ua=[0;0;0]; 
  for m=1:nne         % velocity & derivatives at 8 corner nodes
    mm=ndfn*(m-1);  mm3=mm+1:mm+ndfn;  ncm=nc(m,:);
    S(:,mm3)= JtbD*Z3_V8c{m,k}*Ti{m};  
    Ua = Ua + S(:,mm+1:mm+ndfn)*Vdof(:,m); % A,B,C,u,v,w
  end
  Ub=Jti*Ua;
  UgS=zeros(3,ndfe);
  for m=1:nne         % velocity & derivatives at 8 corner nodes
    mm=ndfn*(m-1);  mm3=mm+1:mm+ndfn; 
    UgS(:,mm3)=JtbD*(Ub(1)*Z3_V8xc{m,k}+Ub(2)*Z3_V8yc{m,k}+Ub(3)*Z3_V8zc{m,k})*Ti{m};
  end 
  Cm = Cm + S'*UgS*W(k)*Det;
  
end    % loop k over quadrature points 

gf=zeros(ndfe,1);
m=0;
for k=1:nne
  m=m+1; gf(m)=nn2nft(1,Elcon(k));  % get global freedom number 
  for k1=2:ndfn
    m=m+1; gf(m)=gf(m-1)+1;  % next 
  end  % if
end %  loop on k
RowNdx=repmat(gf,1,ndfe);
ColNdx=RowNdx';
RowNdx=reshape(RowNdx,ndfe*ndfe,1);
ColNdx=reshape(ColNdx,ndfe*ndfe,1); 
Cm=reshape(Cm,ndfe*ndfe,1);
 
return;

% -----------------------------------------------------------------------

function G=GTL(ni,q,r,s)
% Transposed gradient (derivatives) of scalar trilinear mapping function. 
% The parameter ni can be a vector of coordinate pairs. 
G=[.125*ni(:,1).*(1+ni(:,2).*r).*(1+ni(:,3).*s), .125*ni(:,2).*(1+ni(:,1).*q).*(1+ni(:,3).*s), ... 
    .125*ni(:,3).*(1+ni(:,1).*q).*(1+ni(:,2).*r)];
return;