function Q = ilu_gmres_with_EBC(mat,rhs,EBC,GMRES,Q0,DropTol)
% ILU_GMRES_WITH_EBCx  Solves matrix equation mat*Q = rhs.
% 
% Solves the matrix equation mat*Q = rhs, optionally constrained 
% by Dirichlet boundary conditions described in diri_list, using 
% Matlab's preconditioned gmres sparse solver.  When Dirichlet 
% boundary conditions are provided, the routine enforces them by 
% reordering to partition out Dirichlet degrees of freedom. 
% usage:  Q = ilu_gmres_with_EBC(mat,rhs,EBC,GMRES,Q0,DropTol)
%    or:  Q = ilu_gmres_with_EBC(mat,rhs,EBC,GMRES,Q0)
%    or:  Q = ilu_gmres_with_EBC(mat,rhs,EBC,GMRES)
%    or:  Q = ilu_gmres_with_EBC(mat,rhs,EBC)
%    or:  Q = ilu_gmres_with_EBC(mat,rhs)
%
%  mat  - matrix of linear system to be solved. 
%  rhs  - right hand side of linear system.
%  EBC.dof, EBC.val - (optional) list of Dirichlet boundary 
%          conditions (constraints). May be empty ([]).
%  GMRES - Structure specifying tolerance, max iterations and
%          restarts. Use [] for default values.
%  Q0   - (optional) initial approximation to solution for restart, 
%          may be empty ([]).
%  DropTol - drop tolerance for luinc preconditioner (default=1e-6).
%
% The solution Q is reported with the original ordering restored. 
%
% If specified and not empty, diri_list should have two columns 
%   total and one row for each diri degree of freedom.  The first 
%   column of each row must contain global index of the degree of 
%   freedom.  The second column contains the actual Dirichlet value. 				 
%   If there are no Dirichlet nodes, omit this parameter if not using
%   the restart capabilities, or supply empty matrix ([]) as diri_list.
%
% Nominal values for GMRES for a large difficult problem might be: 
%  GMRES.Tolerance  = 1.e-12, 
%  GMRES.MaxIterates = 75, 
%  GMRES.MaxRestarts = 14.
%
% Jonas Holdeman,   revised February, 2009. 

% Drop tolerance for luinc preconditioner, nominal value - 1.e-6
if (nargin<6 | isempty(DropTol) | DropTol<=0) 
  droptol = 1.e-6;        % default
else droptol = DropTol;   % assigned
end

if nargin<=3 | isempty(GMRES)
   Tol = 1.e-12;  MaxIter = 75;  MaxRstrt = 14;
else
% Tol=tolerance for residual, increase it if solution takes too long.
   Tol=GMRES.Tolerance;
% MaxIter = maximum number of iterations before restart (MT used 5)
   MaxIter=GMRES.MaxIterates;
% MaxRstrt = maximum number of restarts before giving up (MT used 10)
   MaxRstrt=GMRES.MaxRestarts;
end

% check arguments for reasonableness
tdof = size(mat,1);	 % total number of degrees of freedom

% good mat is a square tdof by tdof matrix
if (size(mat,2)~=size(mat,1) | size(size(mat),2)~=2)
  error('mat must be a square matrix')
end

% valid rhs has the dimensions [tdof, 1]
if (size(rhs,1)~=tdof | size(rhs,2)~=1)
  error('rhs must be a column matrix with the same number of rows as mat')
end

% valid dimensions for optional diri_list
if nargin<=2
   EBC=[];
elseif ~isempty(EBC) & (size(EBC.val,1)>=tdof ...
      | size(EBC.dof,1)>=tdof)
   error('check dimensions of EBC')
end
% (optional) valid Q0 is empty or has the dimensions [tdof, 1]
if nargin<=4
   Q0=[];
elseif ~isempty(Q0) & (size(Q0,1)~=tdof | size(rhs,2)~=1)
  error('Q0 must be a column matrix with the same number of rows as mat')
end

% handle the case of no Dirichlet dofs separately
if isempty(EBC)
% skip diri partitioning, solve the system
  % incomplete LU
  [L,U] = lu(mat,droptol);
  Q = gmres(mat,rhs,MaxIter,Tol,MaxRstrt,L,U,Q0);	% GMRES
  
else
% Form list of all DOFs   
  p_vec = [1:tdof]';
% partion out diri dofs
  EBCdofs = EBC.dof(:,1);	 % list of dofs which are Dirichlet
  EBCvals = EBC.val(:,1);  % Dirichlet dof values
  
% form a list of non-diri dofs
  ndro = p_vec(~ismember(p_vec, EBCdofs));	% list of non-diri dofs
  
% Move Dirichlet DOFs to right side
  rhs_reduced = rhs(ndro) - mat(ndro, EBCdofs) * EBCvals;
  
% solve the reduced system (preconditioned gmres)
  A = mat(ndro,ndro);
   
% Compute incomplete LU preconditioner
  [L,U] = lu(A,droptol);      % incomplete LU
   
% Remove Dirichlet DOFs from initial estimate
  if ~isempty(Q0)  Q0=Q0(ndro);  end
 
% solve the reduced system (preconditioned gmres)
  Q_reduced = gmres(A,rhs_reduced,MaxIter,Tol,MaxRstrt,L,U,Q0);
  
% insert the Dirichlet values into the solution
  Q = [Q_reduced; EBCvals];

% calculate p_vec_undo to restore Q to the original dof ordering
  p_vec_undo = zeros(1,tdof);
  p_vec_undo([ndro;EBCdofs]) = [1:tdof];
  
% restore the original ordering
  Q = Q(p_vec_undo);
end