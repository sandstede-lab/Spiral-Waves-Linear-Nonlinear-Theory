function plot_spiral(uout,par,numPar)
% Plots all components

th = linspace(0,2*pi,numPar.nx+1);
r = linspace(par.r1,par.r2,numPar.ny);
[X,Y] = meshgrid(th,r);

figure
for j = 1:par.numVars
    uu = reshape(uout((j-1)*numPar.nx*numPar.ny + 1:j*numPar.nx*numPar.ny),numPar.nx,numPar.ny)';
    u_new = [uu, uu(:,1)]; % Add another angular component using periodic boundary conditions to make sure the full pattern plots
	subplot(1,par.numVars,j)
    [th, r, z] = pol2cart(X,Y,u_new);
    pcolor(th,r,z); shading interp;  % U equation
    title(['U' num2str(j)]);
    colorbar; 
    axis square;
end
drawnow;