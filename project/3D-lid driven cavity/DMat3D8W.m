function [Dm,RowNdx,ColNdx]=DMat3D8W(Xe, Elcon, nn2nft)
% DMAT3D8W - Affine hexahedron (linear 3D Hermite) element diffusion matrix.
%
% 3D linear-complete, normal-conforming, divergence-free, Hermite basis 
%   functions on 8-node rectangular hexahedral elements with 6 DOF per node 
%   using Gauss quadrature on the 2x2x2 reference cube. 
% The assumed nodal numbering starts with 1 at the lower left corner (-1,-1,-1) 
%   of the element . 
%
% Usage:
%   [Dm,Rndx,Cndx] = DMat3D8W(Xe,Elcon,nn2nft)
%   Xe(1,:) -  x-coordinates of 8 corner nodes of element.  
%   Xe(2,:) -  y-coordinates of 8 corner nodes of element.  
%   Xe(3,:) -  z-coordinates of 8 corner nodes of element.  
%   Elcon(8) - connectivity matrix for this element, list of nodes. 
%   nn2nft(1,n) - global freedom number for node n.
%   nn2nft(2,n) - global freedom type for node n.
%
% Calls:
%   V8xyzcW(n,x,y,z).
%
% Jonas Holdeman, July 2011


% Constants and fixed data
nc=[-1,-1,-1; 1,-1,-1; -1,1,-1; 1,1,-1; -1,-1,1; 1,-1,1; -1,1,1; 1,1,1]; % defines corner nodal order
ndfn=6;        % number of degrees of freedom per node. 
nne=8;         % number of nodes per element
ndfe=ndfn*nne;  % number of degrees of freedom per element (48)

% Define 3-point quadrature data once, on first call. 
% Gaussian weights and absissas to integrate 5th degree polynomials exactly on cube. 
global GQC3;
if (isempty(GQC3))       % Define 3-point quadrature data once, on first call. 
   Aq=[-.774596669241483, .000000000000000,.774596669241483]; %Abs
   Hq=[ .555555555555556, .888888888888889,.555555555555556]; %Wts
   GQC3.xa=zeros(27,1); GQC3.ya=zeros(27,1); GQC3.za=zeros(27,1); GQC3.wt=zeros(27,1); 
   nr=0; 
   for nz=1:3; for ny=1:3; for nx=1:3
       nr=nr+1; GQC3.xa(nr)=Aq(nx); GQC3.ya(nr)=Aq(ny); GQC3.za(nr)=Aq(nz); 
       GQC3.wt(nr)=Hq(nx)*Hq(ny)*Hq(nz);
   end; end; end
  GQC3.size=nr; 
end
% Define 4-point quadrature data once, on first call. 
% Gaussian weights and absissas to integrate 7th degree polynomials exactly on cube. 
global GQC4;
if (isempty(GQC4))       % Define 4-point quadrature data once, on first call. 
   Aq=[-.861136311594053,-.339981043584856,.339981043584856, .861136311594053]; %Abs
   Hq=[ .347854845137454, .652145154862546,.652145154862546, .347854845137454]; %Wts
   GQC4.xa=zeros(64,1); GQC4.ya=zeros(64,1); GQC4.za=zeros(64,1); GQC4.wt=zeros(64,1); 
   nr=0; 
   for nz=1:4; for ny=1:4; for nx=1:4
       nr=nr+1; GQC4.xa(nr)=Aq(nx); GQC4.ya(nr)=Aq(ny); GQC4.za(nr)=Aq(nz); 
       GQC4.wt(nr)=Hq(nx)*Hq(ny)*Hq(nz);
   end; end; end
  GQC4.size=nr; 
end
%--------------------------------------------------------------------
%xa=GQC3.xa; ya=GQC3.ya; za=GQC3.za; W=GQC3.wt; Nq=GQC3.size;
xa=GQC4.xa; ya=GQC4.ya; za=GQC4.za; W=GQC4.wt; Nq=GQC4.size;

% ---------------------------------------------------
global Z3_V8xd; global Z3_V8yd; global Z3_V8zd;
if (isempty(Z3_V8xd) | size(Z3_V8xd,2)~=Nq)
  Z3_V8xd=cell(nne,Nq); Z3_V8yd=cell(nne,Nq); Z3_V8zd=cell(nne,Nq);
  for k=1:Nq
    for m=1:nne
      ncm=nc(m,:);
      [Z3_V8xd{m,k},Z3_V8yd{m,k},Z3_V8zd{m,k}]=V8xyzcW(ncm,xa(k),ya(k),za(k));
    end
  end
end   % if (isempty(*))
% ----------------- End fixed data ------------------
%
Ti=cell(nne);
for m=1:nne
  Jt=Xe*GTL(nc(:,:),nc(m,1),nc(m,2),nc(m,3));
  Det=det(Jt);
  JtiD=inv(Jt)*Det;
  J=Jt';
  Ti{m}=blkdiag(J,JtiD);
end  % loop m

%Det=Jt(1,1)*Jt(2,2)*Jt(3,3)-Jt(1,1)*Jt(2,3)*Jt(3,2)-Jt(2,2)*Jt(1,3)*Jt(3,1)...
%   -Jt(3,3)*Jt(1,2)*Jt(2,1)+Jt(1,2)*Jt(2,3)*Jt(3,1)+Jt(2,1)*Jt(1,3)*Jt(3,2);

Dm=zeros(ndfe,ndfe);        % Preallocate arrays
dF1=zeros(3,ndfe); dF2=zeros(3,ndfe); dF3=zeros(3,ndfe);  % derivatives of shape functions

for k=1:Nq  
  Jt=Xe*GTL(nc(:,:),xa(k),ya(k),za(k));    % transpose of Jacobian at (xa,ya,za)
  Det=det(Jt);
%  Det=Jt(1,1)*Jt(2,2)*Jt(3,3)-Jt(1,1)*Jt(2,3)*Jt(3,2)-Jt(2,2)*Jt(1,3)*Jt(3,1)...
%   -Jt(3,3)*Jt(1,2)*Jt(2,1)+Jt(1,2)*Jt(2,3)*Jt(3,1)+Jt(2,1)*Jt(1,3)*Jt(3,2);
  JtbD=Jt/Det;
  Jti=inv(Jt); 
  JtiD=Jti*Det; 
  Ji=Jti';
  for m=1:nne
    mm=6*(m-1);  mm3=mm+1:mm+6;  ncm=nc(m,:);
    dF1(:,mm3)=JtbD*(Ji(1,1)*Z3_V8xd{m,k}+Ji(1,2)*Z3_V8yd{m,k}+Ji(1,3)*Z3_V8zd{m,k})*Ti{m};
    dF2(:,mm3)=JtbD*(Ji(2,1)*Z3_V8xd{m,k}+Ji(2,2)*Z3_V8yd{m,k}+Ji(2,3)*Z3_V8zd{m,k})*Ti{m};
    dF3(:,mm3)=JtbD*(Ji(3,1)*Z3_V8xd{m,k}+Ji(3,2)*Z3_V8yd{m,k}+Ji(3,3)*Z3_V8zd{m,k})*Ti{m};
  end   % loop m  
  
  Dm = Dm + (dF1'*dF1 + dF2'*dF2 + dF3'*dF3)*W(k)*Det;
end  % loop k 

gf=zeros(ndfe,1);
m=0;
for k=1:8
  m=m+1; gf(m)=nn2nft(1,Elcon(k));  % get global freedom number 
  for k1=2:6
    m=m+1; gf(m)=gf(m-1)+1;  % next 
  end  % if
end %  loop on k
RowNdx=repmat(gf,1,ndfe);
ColNdx=RowNdx';
RowNdx=reshape(RowNdx,ndfe*ndfe,1);
ColNdx=reshape(ColNdx,ndfe*ndfe,1); 
Dm=reshape(Dm,ndfe*ndfe,1);
return;

% ---------------------------------------------------------

function G=GTL(ni,q,r,s)
% Transposed gradient (derivatives) of scalar trilinear mapping function. 
% The parameter ni can be a vector of coordinate pairs. 
G=[.125*ni(:,1).*(1+ni(:,2).*r).*(1+ni(:,3).*s), .125*ni(:,2).*(1+ni(:,1).*q).*(1+ni(:,3).*s), ... 
    .125*ni(:,3).*(1+ni(:,1).*q).*(1+ni(:,2).*r)];
return;