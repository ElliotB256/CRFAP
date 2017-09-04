%%
% Check that the coordinate transformations make sense!

RF = 3; % MHz
amp = 0.5 / 0.7; % Gauss
ap = AP.Calculator().LinearPolarised(RF, amp).DontUseParallel();
ap.Atom.F = 1;
disp(ap);

%%
% Check that coupling strength vanishes when parallel to field

sampler = AP.Sampler.LineSampler(ap, pi/2, 0);
sampler.StartB = 2:0.2:4;
sampler.MeshIterations = 3;
sampler.Verbose = 1; sampler.QuadGrad = 100;
sampler.Sort = 0;
sampler.Sample();

plot(sampler.X, sampler.E, '.-');
xlabel('X (\mum)'); ylabel('Eigenenergy (MHz)'); set(gcf, 'Color', 'w');

%%
% Check that coupling is as before on z-axis

sampler = AP.Sampler.LineSampler(ap, 0, 0);
sampler.StartB = 2:0.2:4;
sampler.MeshIterations = 3;
sampler.Verbose = 1; sampler.QuadGrad = 100;
sampler.Sort = 0;
sampler.Sample();

plot(sampler.Z, sampler.E, '.-');
xlabel('Z (\mum)'); ylabel('Eigenenergy (MHz)'); set(gcf, 'Color', 'w');

%%
% Check that coupling is as before on y-axis

sampler = AP.Sampler.LineSampler(ap, 0, pi/2);
sampler.StartB = 2:0.2:4;
sampler.MeshIterations = 3;
sampler.Verbose = 1; sampler.QuadGrad = 100;
sampler.Sort = 0;
sampler.Sample();

plot(sampler.Y, sampler.E, '.-');
xlabel('Y (\mum)'); ylabel('Eigenenergy (MHz)'); set(gcf, 'Color', 'w');