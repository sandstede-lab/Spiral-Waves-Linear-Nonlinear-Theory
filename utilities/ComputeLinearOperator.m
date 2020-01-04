function [L1, L2, D1r] = ComputeLinearOperator(par,numPar)

%% rename parameters
nx = numPar.nx;
ny = numPar.ny;
order = numPar.order;
r1 = par.r1;
r2 = par.r2;

%% theta (angular) direction
switch numPar.thgrid
case 'F' % Fourier
    [~,D2x] = fourdif(nx,2);    % 2nd derivative matrix        
    [~,D1x] = fourdif(nx,1);    % 1st derivative matrix
            
case 'FD' % Finite Differences
    switch order 
    case '2'
    	hx = 2*pi/nx;  % Periodic 
        D1x = sparse(1:nx-1,[2:nx-1 nx],ones(nx-1,1),nx,nx); 
        D1x = D1x - D1x';
        D1x(1,end) = -1; D1x(end,1) = 1; % Periodic boundary conditions
        D1x = D1x/(2*hx);  % 2nd order FD centered difference first derivative matrix

        ex = ones(nx,1);
        D2x = sparse(1:nx-1,[2:nx-1 nx],ones(nx-1,1),nx,nx)- sparse(1:nx,1:nx,ex,nx,nx);
        D2x = D2x + D2x';
        D2x(end,1) = 1;
        D2x(1,end) = 1;
        D2x = D2x/(hx^2);  % 2nd order FD centered difference 2nd derivative matrix

    case '4' % 4th order 
        hx = 2*pi/nx;  % Periodic
        D1x = sparse(1:Nrh-1,[2:Nrh-1 Nrh],8*ones(Nrh-1,1),Nrh,Nrh) - sparse(1:Nrh-2,[3:Nrh-1 Nrh],ones(Nrh-2,1),Nrh,Nrh);
        D1x = (D1x - D1x');
        D1x(1,end-1:end) = [1, -8];D1x(2,end) = 1; D1x(end-1,1) = -1; D1x(end,1:2) = [8,-1]; % Periodic boundary conditions
        D1x = D1x/(12*hx); % First derivative matrix

        D2x = sparse(1:Nrh-1,[2:Nrh-1 Nrh],16*ones(Nrh-1,1),Nrh,Nrh) - sparse(1:Nrh-2,[3:Nrh-1 Nrh],ones(Nrh-2,1),Nrh,Nrh);
        D2x = (D2x + D2x' - 30*speye(Nrh)); 
        D2x(1,end-1:end) = [-1, 16]; D2x(2,end) = -1; D2x(end-1,1) = -1; D2x(end,1:2) = [16,-1]; % Periodic boundary conditions
        D2x = D2x/(12*hx^2); % 2nd derivative matrix
    end
end

%% radial direction
switch numPar.rgrid
case 'FD' 
    [R,D2y,D1y] = Compute_2D_radial_Laplacian_finite_difference(ny,r1,r2,order);
case 'FD_weighted' 
    [R,D2y,~] = Compute_2D_radial_Laplacian_finite_difference_weighted(ny,r2,order);
end

%% identity matrices
Ix  = speye(nx,nx);
Iy = speye(ny,ny);

%% First derivative: d/d(theta)
L1 = kron(Iy,D1x);  % Works for all the cases. Just contains theta derivatives. 

%% Assemble the Laplacian:
L2 = kron(R.^2,D2x) + kron(D2y,Ix);	% Laplace in r
if r1<1e-08		% Add boundary conditions at r=0
	hy = (r2-r1)/(ny-1);
	L2(1:nx,:) = 0;
	L2(1:nx,nx+1:2*nx) = (4/nx)/(hy^2);
	L2(1:nx,1:nx) = -4/(hy^2).*speye(nx);
end
L2 = r2^2*L2;	% Laplace in s=r/r2

%% First derivative in radius if requested
if nargout > 2
   D1r = kron(D1y,Ix); % First derivative matrix for radius    
end
