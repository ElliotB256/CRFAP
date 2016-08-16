%% At the Bottom of the Shell
% Calculates the shape of the potential at the bottom of the shell. What
% happens as we approach the limit that $\omega_r \to 0$?

RF = 4.2;
BRF = 0.6;
BGrad = 400; %G/cm

trap = @(x,z) ShellTrap( x, zeros(size(x)), z, Constants.zeemansplit, BGrad, RF, BRF) + gpe(z, 87);

% Create numerical grid to evaluate the potential over
rew = resonantEllipsoidWidth(RF, BGrad, Constants.zeemansplit);
zr = -rew/2 + (-0.2:0.01:0.5) * rew;
xr = (0:0.002:1.2) * rew;
[x,z] = ndgrid(xr, zr);
tris = delaunay(x,z);
t = trap(x,z);
p = trisurf(tris, x, z, zeros(size(x)), log(t));
shading interp;
axis equal;
xlabel('X ($\mu$m)','Interpreter', 'Latex');
xlabel('Z ($\mu$m)','Interpreter', 'Latex');
view(2)

% Show the location of the minimum
[~,i] = min(trap(x(:),z(:)));
hold on; plot(x(i), z(i), '.', 'MarkerSize', 10, 'Color', 'Green');