%%
trap = @(x,b,z,t) ShellTrap( x, zeros(size(x)), z, zsf, BGrad, 1.4, 0.5);


%%
pot = trap(x,z);
imagesc(pot');

%%
b = 1.5*resonantEllipsoidWidth(1.4, BGrad, zsf);
[x,z] = ndgrid(-b:b/100:b, -b:b/100:b);
tris = delaunay(x,z);
pot = timeAverage(@(t) trap(x,zeros(size(x)),z,t), 10);
p = trisurf(tris, x, z, zeros(size(x)), pot);
shading interp;
axis equal;
xlabel('X ($\mu$m)','Interpreter', 'Latex');
xlabel('Z ($\mu$m)','Interpreter', 'Latex');
view(2)