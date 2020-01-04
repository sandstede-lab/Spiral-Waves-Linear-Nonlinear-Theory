function [D1x,D2x] = ComputeLinearOperator_shortGrid(par,numPar)

%% rename parameters
	nx = numPar.nx;
	ny = numPar.ny;
    order = numPar.order;
	r1 = par.r1;
	r2 = par.r2;
    
% theta direction

switch numPar.thgrid
    case 'F' % Fourier
        [~,D2x] = fourdif(nx,2);    % 2nd derivative matrix        
        [~,D1x] = fourdif(nx,1);    % 1st derivative matrix
            
    case 'FD' % Finite Differences
        switch order 
            case '2'
            hx = 2*pi/nx;  % Periodic so 0 = 2*pi
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

            case '4' % 4th order - Not sure if this is right!
            hx = 2*pi/nx;  % Periodic so 0 = 2*pi
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
Ix  = speye(nx,nx);

% radial direction
switch numPar.rgrid
    case 'FD' % No inner hole - only one that makes sense
        [R,D2y,~] = Compute_2D_radial_Laplacian_finite_difference_nobc_origin(ny,r2,order);
        hy = r2/(ny-1);
        
%     case 'FD_hole' % Inner hole present
%         [r,D2y,~] = Compute_2D_radial_Laplacian_finite_difference_non_zero_radius(ny,r1,r2,order);
%         R = sparse(1:ny,1:ny,1./r,ny,ny);
%     
%     case 'Cheby'  % Chebyshev grid 
%         [r,D2y,~] = Compute_2D_radial_Laplacian_cheb_non_zero_radius(ny,r1,r2);
%         R = sparse(1:ny,1:ny,1./r,ny,ny);
%         
%     case 'FD_symm' % No inner hole, assumes radial symmetry (don't use for spiral)
%         [R,D2y,~] = Compute_2D_radial_Laplacian_finite_difference_symm(ny,r2);
        
	
    
end
Iy = speye(ny,ny);

% First derivative: d/d(theta)
D1x = kron(Iy,D1x);  % Works for all the cases. Just contains theta derivatives. 
D1x = D1x(nx:end,nx:end); % Remove r=0 except for 1 row and comlumn
D1x(:,1) = 0; D1x(1,:) = 0;

% Assemble the Laplacian: 
L1 = kron(R.^2,D2x);    % Components of the Laplacian
L1 = L1(nx:end,nx:end); % r1, r2 do not need to be updated
L1(:,1) = 0; L1(1,:) = 0; % r0 component

L2 = kron(D2y,Ix);      % r1, r2 blocks need to be updated. 
L2 = L2(nx:end,nx:end);
L2(:,1) = 0; L2(1,:) = 0;   %r0 component

% Radius 1 update
L2(2:nx+1, 2:nx+1) =  -(5/2).*(1./hy.^2).*speye(nx);    % from d_rr. Note that the off-diagonal terms of (1/r)*d_r and d_rr cancel!
L2(2:nx+1,1) = 2./(3.*hy.^2);
L2(nx+2:2*nx+1, 1) = -1./(24.*hy.^2);

% Radius 2 update  - don't need to do anythin! Terms cancel or are covered
% by Kron!

% 9 point laplacian at the origin
%D2x  =  kron(R.^2,D2x) + kron(D2y,Ix); 
D2x = L1 + L2;
tmp = zeros(1,2*nx);
tmp([1,floor(nx/4)+1, floor(nx/2)+1, floor(3*nx/4) + 1]) = 4/3;        % radius 1 - best if nx is multiple of 4
tmp([nx+1,floor(5/4*nx)+1 floor(3*nx/2)+1, floor(7*nx/4)+1]) = -1/12;  % radius 2
D2x(1,1:2*nx+1) = [-5, tmp]./(hy^2);










