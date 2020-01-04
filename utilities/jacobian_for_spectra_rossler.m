function [f,J] = jacobian_for_spectra_rossler(u,L1,L2,par,numPar)
% J = Jacobian on the reduced polar grid, with only one grid point at the
% origin
% Need u, L1, L2 in short grid format

nx = numPar.nx;
ny = numPar.ny;

a = par.a;
b = par.b;
c = par.c;
omega = par.omega;

delta1 = par.delta1;
delta2 = par.delta2;
delta3 = par.delta3;

% Variables: Put into short grid format - only one point at the origin.
U = u(1:nx*ny); U = U(nx:end);
V = u(nx*ny+1:2*nx*ny); V = V(nx:end);
W = u(2*nx*ny+1:3*nx*ny); W = W(nx:end);


% Reaction terms
fU = -V - W;
fV = U + a*V;
fW = U.*W - c.*W + b;

line1 = delta1.*(L2*U) + omega.*(L1*U) + fU;
line2 = delta2.*(L2*V) + omega.*(L1*V) + fV;
line3 = delta3.*(L2*W) + omega.*(L1*W) + fW;

f = [line1; line2; line3];  % Not really used, just for completeness
   
% Jacobian              
if (nargout > 1)
    
    fU_V = -1;  fU_W = -1;
    fV_U = 1;   fV_V = a;
    fW_U = W;   fW_W = U - c;
    
    I = speye(nx*(ny-1)+1,nx*(ny-1)+1);
    
    
    dU = [delta1.*L2 + omega.*L1;
        fV_U.*I;
        spdiags(fW_U,0,nx*(ny-1)+1,nx*(ny-1)+1)];
    
    dV = [fU_V.*I;
        delta2.*L2 + omega.*L1 + fV_V.*I;
        sparse(nx*(ny-1)+1,nx*(ny-1)+1)];
    
    dW = [fU_W.*I;
        sparse(nx*(ny-1)+1,nx*(ny-1)+1);
        delta3.*L2 + omega.*L1 + spdiags(fW_W,0,nx*(ny-1)+1,nx*(ny-1)+1)];
    
        
    J = [dU, dV, dW];
    J = sparse(J);
   
          
end




