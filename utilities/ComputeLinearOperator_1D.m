function [L1,L2] = ComputeLinearOperator_1D(par,numPar)
% Only does periodic boundary conditions
%% rename parameters
	nx = numPar.nx;
    order = numPar.order;
    Lx = numPar.Lx;


switch numPar.thgrid
    case 'F' % Fourier
        [~,L2] = fourdif(nx,2);    % 2nd derivative matrix        
        [~,L1] = fourdif(nx,1);    % 1st derivative matrix
            
    case 'FD' % Finite Differences
        switch order 
            case '2'
            hx = Lx/nx;  % Periodic so 0 = Lx
            L1 = sparse(1:nx-1,[2:nx-1 nx],ones(nx-1,1),nx,nx); 
            L1 = L1 - L1';
            L1(1,end) = -1; L1(end,1) = 1; % Periodic boundary conditions
            L1 = L1/(2*hx);  % 2nd order FD centered difference first derivative matrix

            ex = ones(nx,1);
            L2 = sparse(1:nx-1,[2:nx-1 nx],ones(nx-1,1),nx,nx)- sparse(1:nx,1:nx,ex,nx,nx);
            L2 = L2 + L2';
            L2(end,1) = 1;
            L2(1,end) = 1;
            L2 = L2/(hx^2);  % 2nd order FD centered difference 2nd derivative matrix

            case '4' % 4th order
            hx = Lx/nx;  % Periodic so 0 = Lx
            L1 = sparse(1:nx-1,[2:nx-1 nx],8*ones(nx-1,1),nx,nx) - sparse(1:nx-2,[3:nx-1 nx],ones(nx-2,1),nx,nx);
            L1 = (L1 - L1'); 
            L1(1,end-1:end) = [1, -8]; L1(2,end) = 1; L1(end-1,1) = -1; L1(end,1:2) = [8,-1]; % Periodic boundary conditions
            L1 = L1/(12*hx); % First derivative matrix

            L2 = sparse(1:nx-1,[2:nx-1 nx],16*ones(nx-1,1),nx,nx) - sparse(1:nx-2,[3:nx-1 nx],ones(nx-2,1),nx,nx);
            L2 = (L2 + L2' - 30*speye(nx)); 
            L2(1,end-1:end) = [-1, 16]; L2(2,end) = -1; L2(end-1,1) = -1; L2(end,1:2) = [16,-1]; % Periodic boundary conditions
            L2 = L2/(12*hx^2); % 2nd derivative matrix
        end
        
    case 'FD_Neumann'   % 4th order Neumann boundary conditions
        
        hx = Lx/(nx-1);
        L1 = sparse(1:nx-1,[2:nx-1 nx],8*ones(nx-1,1),nx,nx) - sparse(1:nx-2,[3:nx-1 nx],ones(nx-2,1),nx,nx);
        L1 = (L1 - L1')/12; L1(1:2,:) = 0; L1(end-1:end,:) = 0; L1(2,1:3) = [-1/2, 0, 1/2]; L1(end-1,end-2:end) = [-1/2, 0, 1/2]; % Use 2nd order at boundary
        L1 = L1/hx; % First derivative matrix

        L2 = sparse(1:nx-1,[2:nx-1 nx],16*ones(nx-1,1),nx,nx) - sparse(1:nx-2,[3:nx-1 nx],ones(nx-2,1),nx,nx);
        L2 = (L2 + L2' - 30*speye(nx))/12; 
        L2(1:2,:) = 0; L2(end-1:end,:) = 0; L2(2,1:3) = [1, -2, 1];  % Neumann boundary conditions: use 2nd order at boundary
        L2(end-1,end-2:end) = [1, -2, 1]; L2(1,1:2) = [-2,2]; L2(end,end-1:end) = [2,-2];
        L2 = L2/hx^2; % 2nd derivative matrix   
        
    case 'FD_Dirichlet_Neumann'   % 4th order Neumann boundary conditions
        
        hx = Lx/(nx-1);
        L1 = sparse(1:nx-1,[2:nx-1 nx],8*ones(nx-1,1),nx,nx) - sparse(1:nx-2,[3:nx-1 nx],ones(nx-2,1),nx,nx);
        L1 = (L1 - L1')/12; L1(1:2,:) = 0; L1(end-1:end,:) = 0; L1(1,1:2) = [0, 1/2]; L1(2,1:3) = [-1/2, 0, 1/2]; L1(end-1,end-2:end) = [-1/2, 0, 1/2]; % Use 2nd order at boundary. Keeping first one as is
        L1 = L1/hx; % First derivative matrix

        L2 = sparse(1:nx-1,[2:nx-1 nx],16*ones(nx-1,1),nx,nx) - sparse(1:nx-2,[3:nx-1 nx],ones(nx-2,1),nx,nx);
        L2 = (L2 + L2' - 30*speye(nx))/12; 
        L2(1:2,:) = 0; L2(end-1:end,:) = 0; L2(2,1:3) = [1, -2, 1];  % Neumann boundary conditions: use 2nd order at boundary. Keeping left endpoint as is
        L2(end-1,end-2:end) = [1, -2, 1]; L2(1,1:2) = [-2,1]; L2(end,end-1:end) = [2,-2];
        L2 = L2/hx^2; % 2nd derivative matrix   
    
        
end



