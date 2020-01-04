function [f,J] = Barkley_weighted_operator(U,V,L1,L2,L1r,par,numPar)
% Weighted opeartor for spiral wave on radial disk 
% Want opeartors in the short grid format

nx = numPar.nx;
ny = numPar.ny;

a = par.a;
b = par.b;
ep = par.ep;
delta = par.delta;  % Diffusion coefficient of v-equation
w = par.w;      % Weight
omega = par.omega;


hxn = par.r2/(ny-1);
R = (1:ny-1)*hxn;                      % radial mesh (not including origin)
R = repmat(R,nx,1);
R = R(:); R = [1;R];            % Add one more point for the origin
R = w./R;

% Barkley model (weighted part in () )
nonlin = (1/ep)*U.*(1-U).*(U - (V+b)/a);
line1 =   L2*U + (  w.^2.*U + 2.*w.*(L1r*U) + R.*U  )  + omega*(L1*U)  + nonlin;
line2 =  delta.*L2*V + delta.*(  w.^2.*V + 2.*w.*(L1r*V)+  R.*V  ) + omega*L1*V + U-V;
f = [line1; line2];     % function...not actually used... just to help take jacobian
   
% Jacobian              
if (nargout > 1)
nonlin_jacob1 = (1/ep)*(2*(1+b/a)*U-b/a-(1/a)*V-3*U.^2+(2/a)*U.*V);
nonlin_jacob2 = (1/ep)*(-(1/a)*U + (1/a)*U.^2);


dU = [ L2  + 2.*w.*L1r + omega*L1  + spdiags(nonlin_jacob1 + R +  w.^2,0,nx*(ny-1)+1,nx*(ny-1)+1);
    speye(nx*(ny-1)+1,nx*(ny-1)+1)];

dV = [spdiags(nonlin_jacob2,0,nx*(ny-1)+1,nx*(ny-1)+1);
    delta.*L2 + 2.*delta.*w.*L1r + omega*L1 + spdiags(-1 + delta.*R + delta.*w.^2,0,nx*(ny-1)+1,nx*(ny-1)+1)];


J = [dU, dV];

% J = [-omega*L1 + L2 + spdiags(nonlin_jacob1,0,nx*ny,nx*ny), spdiags(nonlin_jacob2,0,nx*ny,nx*ny), -L1*U;
%           spdiags(ones(nx*ny,1),0,nx*ny,nx*ny), spdiags(-ones(nx*ny,1),0,nx*ny,nx*ny) - omega*L1 + delta*L2, -L1*V;
%           phase_jacob, sparse(1,nx*ny), 0];
   
          
end


    
    
    
    
    
    
    

