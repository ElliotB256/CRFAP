%% Test Coupling Strength
% Renders numerous dressed atom scenarios as a test that the code is
% functioning correctly. Because the coupling strength equals the
% difference between two adjacent dressed energy levels in the same
% manifold, we can quite straightforwardly measure it from the calculate
% eigenenergies.

%% Test 1: Linear, X
% We expect to see coupling strength minimum on the intersection with the
% x-axis

RF = 3; % MHz
amp = 0.5 / 0.7; % Gauss
ap = AP.Calculator().LinearPolarised(RF, amp).DontUseParallel();
ap.Atom.F = 1;
ap.BX = amp; ap.BY = 0; ap.BZ = 0;

sampler = AP.Sampler.QuadrupoleSurfaceSampler(ap);
sampler.QuadGrad = 100;
sampler.Verbose = 1;
sampler.Theta = linspace(0, pi, 20);
sampler.Gamma = linspace(0, 2*pi, 20);
sampler.Sample();

E = sampler.Eigenenergies;
CS = E(2,:) - E(1,:);

clf; set(gcf, 'Color', 'w');
th = sampler.GetThetas();
ga = sampler.GetGammas();

% grid data, so that we can contour it
[gaG,thG] = meshgrid(sampler.Gamma, sampler.Theta);
csG = griddata(ga,th,CS,gaG,thG);

tris = delaunay(ga, th);
trisurf(tris, ga, th, zeros(size(CS)), CS); shading interp
hold on; contour(gaG, thG, csG, 'k'); hold off;
view(0,90); axis tight;
xlim([ 0 2*pi]); set(gca, 'XTick', [ 0 pi 2*pi ], 'XTickLabel', { '0', '\pi', '2 \pi' });
ylim([ 0 1*pi]); set(gca, 'YTick', [ 0 pi ], 'YTickLabel', { '0', '\pi' });

%% Test 2: Linear, Y
% We expect to see coupling strength minimum on the intersection with the
% y-axis

RF = 3; % MHz
amp = 0.5 / 0.7; % Gauss
ap = AP.Calculator().LinearPolarised(RF, amp).DontUseParallel();
ap.Atom.F = 1;
ap.BX = 0; ap.BY = amp; ap.BZ = 0;

sampler = AP.Sampler.QuadrupoleSurfaceSampler(ap);
sampler.QuadGrad = 100;
sampler.Verbose = 1;
sampler.Theta = linspace(0, pi, 20);
sampler.Gamma = linspace(0, 2*pi, 20);
sampler.Sample();

E = sampler.Eigenenergies;
CS = E(2,:) - E(1,:);

clf; set(gcf, 'Color', 'w');
th = sampler.GetThetas();
ga = sampler.GetGammas();

% grid data, so that we can contour it
[gaG,thG] = meshgrid(sampler.Gamma, sampler.Theta);
csG = griddata(ga,th,CS,gaG,thG);

tris = delaunay(ga, th);
trisurf(tris, ga, th, zeros(size(CS)), CS); shading interp
hold on; contour(gaG, thG, csG, 'k'); hold off;
view(0,90); axis tight;
xlim([ 0 2*pi]); set(gca, 'XTick', [ 0 pi 2*pi ], 'XTickLabel', { '0', '\pi', '2 \pi' });
ylim([ 0 1*pi]); set(gca, 'YTick', [ 0 pi ], 'YTickLabel', { '0', '\pi' });

%% Test 3: Linear, Z
% We expect to see coupling strength minimum on the intersection with the
% z-axis

RF = 3; % MHz
amp = 0.5 / 0.7; % Gauss
ap = AP.Calculator().LinearPolarised(RF, amp).DontUseParallel();
ap.Atom.F = 1;
ap.BX = 0; ap.BY = 0; ap.BZ = amp;

sampler = AP.Sampler.QuadrupoleSurfaceSampler(ap);
sampler.QuadGrad = 100;
sampler.Verbose = 1;
sampler.Theta = linspace(0, pi, 20);
sampler.Gamma = linspace(0, 2*pi, 20);
sampler.Sample();

E = sampler.Eigenenergies;
CS = E(2,:) - E(1,:);

clf; set(gcf, 'Color', 'w');
th = sampler.GetThetas();
ga = sampler.GetGammas();

% grid data, so that we can contour it
[gaG,thG] = meshgrid(sampler.Gamma, sampler.Theta);
csG = griddata(ga,th,CS,gaG,thG);

tris = delaunay(ga, th);
trisurf(tris, ga, th, zeros(size(CS)), CS); shading interp
hold on; contour(gaG, thG, csG, 'k'); hold off;
view(0,90); axis tight;
xlim([ 0 2*pi]); set(gca, 'XTick', [ 0 pi 2*pi ], 'XTickLabel', { '0', '\pi', '2 \pi' });
ylim([ 0 1*pi]); set(gca, 'YTick', [ 0 pi ], 'YTickLabel', { '0', '\pi' });

%% Test 4: Circular
% We expect to see coupling strength minimum at the top of the spheroid

RF = 3; % MHz
amp = 0.5 / 0.7; % Gauss
ap = AP.Calculator().CircularPolarised(RF, amp).DontUseParallel();
ap.Atom.F = 1;

sampler = AP.Sampler.QuadrupoleSurfaceSampler(ap);
sampler.QuadGrad = 100;
sampler.Verbose = 1;
sampler.Theta = linspace(0, pi, 20);
sampler.Gamma = linspace(0, 2*pi, 20);
sampler.Sample();

E = sampler.Eigenenergies;
CS = E(2,:) - E(1,:);

clf; set(gcf, 'Color', 'w');
th = sampler.GetThetas();
ga = sampler.GetGammas();

% grid data, so that we can contour it
[gaG,thG] = meshgrid(sampler.Gamma, sampler.Theta);
csG = griddata(ga,th,CS,gaG,thG);

tris = delaunay(ga, th);
trisurf(tris, ga, th, zeros(size(CS)), CS); shading interp
hold on; contour(gaG, thG, csG, 'k'); hold off;
view(0,90); axis tight;
xlim([ 0 2*pi]); set(gca, 'XTick', [ 0 pi 2*pi ], 'XTickLabel', { '0', '\pi', '2 \pi' });
ylim([ 0 1*pi]); set(gca, 'YTick', [ 0 pi ], 'YTickLabel', { '0', '\pi' });

%% Test 4: Circular, F=2
% We expect to see coupling strength minimum at the top of the spheroid

RF = 3; % MHz
amp = 0.5 / 0.7; % Gauss
ap = AP.Calculator().CircularPolarised(RF, amp).DontUseParallel().OfSpecies('Species', 87, 'F', 2);
ap.Atom.F = 1;

sampler = AP.Sampler.QuadrupoleSurfaceSampler(ap);
sampler.QuadGrad = 100;
sampler.Verbose = 1;
sampler.Theta = linspace(0, pi, 20);
sampler.Gamma = linspace(0, 2*pi, 20);
sampler.Sample();

E = sampler.Eigenenergies;
CS = E(2,:) - E(1,:);

clf; set(gcf, 'Color', 'w');
th = sampler.GetThetas();
ga = sampler.GetGammas();

% grid data, so that we can contour it
[gaG,thG] = meshgrid(sampler.Gamma, sampler.Theta);
csG = griddata(ga,th,CS,gaG,thG);

tris = delaunay(ga, th);
trisurf(tris, ga, th, zeros(size(CS)), CS); shading interp
hold on; contour(gaG, thG, csG, 'k'); hold off;
view(0,90); axis tight;
xlim([ 0 2*pi]); set(gca, 'XTick', [ 0 pi 2*pi ], 'XTickLabel', { '0', '\pi', '2 \pi' });
ylim([ 0 1*pi]); set(gca, 'YTick', [ 0 pi ], 'YTickLabel', { '0', '\pi' });

%%
% Plot a resonant spheroid with coupling strength as color

clf; set(gcf, 'Color', 'w');

x = sampler.X;
y = sampler.Y;
z = sampler.Z;

tris = delaunay(x,y,z);
trisurf(tris,x,y,z,CS); shading interp
xlabel('X (\mum)'); ylabel('Y (\mum)'); zlabel('Z (\mum)');