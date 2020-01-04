function [f,J] = Barkley_weighted_operator(u, L1, L2, L1r, par, numPar)
% Weighted operator for spiral wave on radial disk 
% Want operators in the short grid format

nx = numPar.nx;
ny = numPar.ny;

a = par.a;
b = par.b;
ep = par.ep;
delta = par.delta;
w = par.w;					% Weight
omega = par.omega;

hxn = par.r2/(ny-1);
R = (1:ny-1)*hxn;			% radial mesh (not including origin)
R = repmat(R,nx,1);
R = R(:); R = [1;R];		% Add one more point for the origin
R = w./R;

% Variables: Put into short grid format - only one point at the origin.
U = u(1:nx*ny); U = U(nx:end);
V = u(nx*ny+1:2*nx*ny); V = V(nx:end);

% Barkley model (weighted part in () )
nl1 = (1/ep)*U.*(1-U).*(U - (V+b)/a);
nl2 = (1-6.75.*U.*(U-1).^2).*(U<=1).*(U>=1/3) + 1.*(U>1);
line1 = L2*U       + omega*L1*U + nl1     + (w.^2.*U + 2.*w.*(L1r*U) + R.*U) ;
line2 = delta*L2*V + omega*L1*V + nl2 - V + delta.*(w.^2.*V + 2.*w.*(L1r*V) + R.*V);
f = [line1; line2];

% Jacobian              
if (nargout > 1)
	I = speye(nx*(ny-1)+1,nx*(ny-1)+1);

	fU_U = (1/ep)*(2*(1+b/a)*U-b/a-(1/a)*V-3*U.^2+(2/a)*U.*V);
	fU_V = (1/ep)*(-(1/a)*U + (1/a)*U.^2);
	fV_U = -6.75.*((U-1).^2 +2.*U.*(U-1)).*(U<=1).*(U>=1/3);
	fV_V = -1;

	dU = [ L2 + 2.*w.*L1r + omega*L1 + spdiags(fU_U + R + w.^2,0,nx*(ny-1)+1,nx*(ny-1)+1);
		fV_U.*I];

	dV = [spdiags(fU_V,0,nx*(ny-1)+1,nx*(ny-1)+1);
		delta.*L2 + 2.*delta.*w.*L1r + omega*L1 + spdiags(fV_V + delta.*R + delta.*w.^2,0,nx*(ny-1)+1,nx*(ny-1)+1)];

	J = [dU, dV];
end
