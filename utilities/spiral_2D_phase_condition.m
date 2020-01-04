
% Compute the linear operator: L1 is d/d(theta), L2 is radial laplacian
[L1, L2] = ComputeLinearOperator(par,numPar); 
numPar.L1 = L1; 
numPar.L2 = L2;

% Set up the phase condition for the spiral 
pc = zeros(numPar.nx,numPar.ny);
pc(:,ceil(numPar.ny/2)) = 1;   % elements corresponding to R/2 circle are 1
phase_cond.pc = pc(:);

us = (phase_cond.pc).*U0(1:numPar.nx*numPar.ny);  % Select R/2 elements from initial condition.
phase_cond.u_star = us(phase_cond.pc>0);

us_th = L1*U0(1:numPar.nx*numPar.ny); % d/d(theta) of phase condition
phase_cond.u_star_th = us_th(phase_cond.pc > 0);
