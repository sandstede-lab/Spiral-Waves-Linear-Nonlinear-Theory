%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute 2D-radial finite-difference Laplacian - Neumann bcs on r=[0,Lrh]
% Inputs: Lrh - length of domain, Nrh - number of mesh points
% Outputs: r - mesh, L - Laplacian, Dx - 1D differentiation matrix,
% Note: returned matrices are sparse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rn,L,Dxn] = Compute_2D_radial_Laplacian_finite_difference_nobc_origin(Nrh,Lrh,order)

hxn = Lrh/(Nrh-1);
rn = (0:Nrh-1)'*hxn; rn(1) = 1;                      % radial mesh
Rn = sparse(1:Nrh,1:Nrh,1./rn,Nrh,Nrh); 

%% Finite-difference matrices
switch order
    case '2'
    
    ex = ones(Nrh,1);

    Dxn = sparse(1:Nrh-1,[2:Nrh-1 Nrh],ones(Nrh-1,1)/2,Nrh,Nrh); % 1st derivative matrix
    Dxn = (Dxn - Dxn')/hxn; 
    Dxn(1,2) = 0; Dxn(Nrh,Nrh-1) = 0; % Neumann boundary conditions

    D2xn = sparse(1:Nrh-1,[2:Nrh-1 Nrh],ones(Nrh-1,1),Nrh,Nrh) - sparse(1:Nrh,1:Nrh,ex,Nrh,Nrh); % 2nd derivative matrix
    D2xn = (D2xn + D2xn'); 
    D2xn(1,1)=2; D2xn(1,2) = -2;
    D2xn(Nrh,Nrh)=-2;     % Neumann boundary conditions - alter r = 0 conditions after kron product
    D2xn(Nrh,Nrh-1) = 2;
    D2xn = D2xn/hxn^2;
    
    case '4'
    Dxn = sparse(1:Nrh-1,[2:Nrh-1 Nrh],8*ones(Nrh-1,1),Nrh,Nrh) - sparse(1:Nrh-2,[3:Nrh-1 Nrh],ones(Nrh-2,1),Nrh,Nrh);
    Dxn = (Dxn - Dxn')/12;  Dxn(end-1:end,:) = 0;  Dxn(end-1,end-2:end) = [-1/2,0,1/2]; % 2nd order at outer boundary. r = 0 will be altered after kron
    Dxn = Dxn/hxn; % First derivative matrix
    

    D2xn = sparse(1:Nrh-1,[2:Nrh-1 Nrh],16*ones(Nrh-1,1),Nrh,Nrh) - sparse(1:Nrh-2,[3:Nrh-1 Nrh],ones(Nrh-2,1),Nrh,Nrh);
    D2xn = (D2xn + D2xn' - 30*speye(Nrh))/12; 
     D2xn(end-1:end,:) = 0;  % Neumann boundary conditions: use 2nd order at outer boundary. r = 0 will be altered after kron
    D2xn(end-1,end-2:end) = [1, -2, 1];  D2xn(end,end-1:end) = [2,-2];
    D2xn = D2xn/hxn^2; % 2nd derivative matrix   
    
end
        

%% 2D radial
L = D2xn + Rn*Dxn; % d_rr + d_r/r


