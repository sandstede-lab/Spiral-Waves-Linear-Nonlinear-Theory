function [f,J]  = Rossler_essSpec_fcn(u,u_infty,L1,L2,par,numPar,phase_cond)
% Spiral dispersion relation with free parameter the temporal eigenvalue.
% Used in computation of essential spectrum
% In spiral coordinates


nx = numPar.nx;
ny = numPar.ny;

a = par.a;
b = par.b;
c = par.c;
delta1 = par.delta1;
delta2 = par.delta2;
delta3 = par.delta3;
kappa = par.kappa;
nu = par.nu;
omega = par.omega;

U = u(1:nx*ny);
V = u(nx*ny+1:2*nx*ny);
W = u(2*nx*ny+1:3*nx*ny);
lambda = u(end); % Free parameter = temporal eigenvalue, lambda


U_infty = u_infty(1:nx); % 1D asymptotic wave train
V_infty = u_infty(nx+1:2*nx);
W_infty = u_infty(2*nx+1:3*nx);

% Phase Condtiion
u_phase = phase_cond.u_old;

% Linearizations of non-linear terms around the wave trains
fU_V = -1;  fU_W = -1;
fV_U = 1;   fV_V = a;
fW_U = W_infty;   fW_W = U_infty - c;
    

% spiral dispersion relation: nu is a COMPLEX number
line1 = kappa^2*delta1*(L2*U) + 2*delta1*kappa*nu*(L1*U) + delta1*nu^2*U + omega*(L1*U) + fU_V.*V + fU_W.*W - lambda*U;
line2 = kappa^2*delta2*(L2*V) + 2*delta2*kappa*nu*(L1*V) + delta2*nu^2*V + omega*(L1*V) + fV_U.*U + fV_V.*V - lambda*V;
line3 = kappa^2*delta3*(L2*W) + 2*delta3*kappa*nu*(L1*W) + delta3*nu^2*W + omega*(L1*W) + fW_U.*U + fW_W.*W - lambda*W;
line4 = (u_phase'*[U;V;W]) - 1; 

f = [line1; line2;line3;line4];

% Jacobian              
if (nargout > 1)
    
    I = speye(nx,nx);
    
    dU1 = kappa^2.*delta1.*L2 + 2.*delta1*kappa*nu*L1  +omega*L1 + (delta1*nu^2 - lambda).*I;
    dU2 = fV_U.*I;
    dU3 = spdiags(fW_U,0,nx,nx);
    dU4 = u_phase(1:nx)';
    
    dU = [dU1; dU2; dU3; dU4];
    
    dV = [fU_V.*I;
        kappa^2*delta2*L2 + 2*delta2*kappa*nu*L1  + omega*L1 + (delta2*nu^2 + fV_V - lambda).*I;
        sparse(nx,nx);
        u_phase(nx+1:2*nx)'];
    
    dW = [fU_W.*I;
        sparse(nx,nx);
        kappa^2*delta3*L2 + 2*delta3*kappa*nu*L1  + omega*L1 + spdiags(delta3*nu^2 + fW_W - lambda,0,nx,nx);
        u_phase(2*nx+1:3*nx)'];

   dlambda = [-U;
        -V;
        -W;
        0];
   
   J = [dU, dV, dW, dlambda];
   
   J = sparse(J);
end



    
    
    
    
    
    
    