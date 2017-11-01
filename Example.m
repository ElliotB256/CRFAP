%% Examples
% A list of examples showing how to use the AP package.

%% Create an Adiabatic Potential (AP) Calculator object
% The APCalculator object configures the rf-dressed atom system. This
% includes describing the dressing fields, and specifying properties of the
% atom. Create the AP object, then either set properties directly or use
% the fluent interfacing syntax (as below) to configure eg. circular
% polarised dressing rf.
% 
% After configuring the system, we can get dressed energies using
% GetDressedEnergies.

RF = 3; % MHz
amp = 0.5; % Gauss
ap = AP.Calculator().CircularPolarised(RF, amp).DontUseParallel();
disp(ap);

zeemanSplitting = 2:0.1:4;
E = ap.GetDressedEnergies(zeemanSplitting);
plot(zeemanSplitting, E, '.-');
xlabel('Zeeman splitting (MHz)'); ylabel('Eigenenergy (MHz)'); set(gcf, 'Color', 'w');

%% 
% ...or use linear polarised RF:

RF = 3; % MHz
amp = 0.5; % Gauss
ap = AP.Calculator().LinearPolarised(RF, amp).DontUseParallel();
disp(ap);

zeemanSplitting = 2:0.1:4;
E = ap.GetDressedEnergies(zeemanSplitting);
plot(zeemanSplitting, E, '.-');
xlabel('Zeeman splitting (MHz)'); ylabel('Eigenenergy (MHz)'); set(gcf, 'Color', 'w');

%% Use of the LineSampler
% In this example we use a LineSampler to calculate dressed eigenenergies
% along the z-axis, rather than specify the points to calculate manually.
% This allows us to take advantage of support for meshing and sorting the
% eigenenergies.

RF = 3; % MHz
amp = 0.5 / 0.7; % Gauss
ap = AP.Calculator().CircularPolarised(RF, amp).DontUseParallel();
ap.Atom.F = 1;
disp(ap);

sampler = AP.Sampler.LineSampler(ap);
sampler.StartB = 2:0.2:4;
sampler.MeshIterations = 5;
sampler.Verbose = 1; sampler.QuadGrad = 100;
sampler.Sample();

plot(sampler.Z, sampler.E, '.-');
xlabel('Z (\mum)'); ylabel('Eigenenergy (MHz)'); set(gcf, 'Color', 'w');

%% MRF: Simple example
% An example of how to implement multiple RF fields.
RF  = [ 3.6 3.8 4.0 ]; % MHz
amp = [ 0.16 0.2 0.16 ] / 0.7; % Gauss
% 
ap = AP.Calculator().LinearPolarised(RF, amp).DontUseParallel().OfSpecies('species', 87, 'F', 1);
ap.PY = pi/2; ap.BY = ap.BX * 0.2;
sampler = AP.Sampler.LineSampler(ap);
sampler.StartB = linspace(3.5, 4.1, 10);
sampler.MeshIterations = 3;
sampler.Verbose = 1;
sampler.Sample();

clf; plot(sampler.B, sampler.E, '.-');
xlabel('Zeeman splitting (MHz)'); ylabel('Eigenenergy (MHz)'); set(gcf, 'Color', 'w');

%% Mixtures of species
% You can use the code to also calculate potentials for different species:

RF = [ 3 4.5 ]; % MHz
amp = [ 0.5 0.5 ]; % Gauss
ap87 = AP.Calculator().CircularPolarised(RF, amp).DontUseParallel().OfSpecies('Species', 87, 'F', 1);
ap85 = AP.Calculator().CircularPolarised(RF, amp).DontUseParallel().OfSpecies('Species', 85, 'F', 2);

sampler87 = AP.Sampler.LineSampler(ap87);
sampler87.StartB = 2:0.2:6;
sampler87.MeshIterations = 5;
sampler87.Verbose = 1; sampler87.QuadGrad = 100;
sampler87.Sample();

sampler85 = AP.Sampler.LineSampler(ap85);
sampler85.StartB = 2:0.2:6;
sampler85.MeshIterations = 5;
sampler85.Verbose = 1; sampler85.QuadGrad = 100;
sampler85.Sample();

plot(sampler85.Z, sampler85.E, '.', 'Color', [  0.0000    0.5156    0.7813 ]); hold on;
plot(sampler87.Z, sampler87.E, '.', 'Color', [  0.3203    0.1758    0.3359 ]); hold off;
xlabel('Z (\mum)'); ylabel('Eigenenergy (MHz)'); set(gcf, 'Color', 'w');

%% MRF: Arbitrary polarisation
% In this more complicated example, we use a multiple RF field to
% independently manipulate F=1 and F=2 atoms. We specify an elliptical RF
% dressed field, pushing the barrier down for F=1 atoms and up for F=2 atoms.

RF  = [ 3.6 3.8 4.0 ]; % MHz
amp = [ 0.16 0.2 0.16 ] / 0.7; % Gauss
 
factor = 0.4;
ap = AP.Calculator().LinearPolarised(RF, amp).DontUseParallel().OfSpecies('species', 87, 'F', 1);
ap.PY = [ 0 pi/2 0 ]; ap.BY = ap.BX .* factor .* [ 0 1 0 ]';

sampler = AP.Sampler.LineSampler(ap);
sampler.StartB = linspace(3.4, 4.2, 15);
sampler.MeshIterations = 3;
sampler.Verbose = 1;
sampler.Sample();

clf; h1 = plot(sampler.B, sampler.E, '.-r'); hold on;


ap = AP.Calculator().LinearPolarised(RF, amp).DontUseParallel().OfSpecies('species', 87, 'F', 2);
ap.PY = [ 0 pi/2 0 ]; ap.BY = ap.BX * factor .* [ 0 1 0 ]';

sampler = AP.Sampler.LineSampler(ap);
sampler.StartB = linspace(3.4, 4.2, 15);
sampler.MeshIterations = 3;
sampler.Verbose = 1;
sampler.Sample();

h2 = plot(sampler.B, sampler.E, '.-k'); hold off;

axis tight;
xlabel('Zeeman splitting (MHz)'); ylabel('Eigenenergy (MHz)'); set(gcf, 'Color', 'w');
legend([h1(1) h2(1)], 'F = 1', 'F = 2');

%% Plane Sampler
% Use the Plane Sampler to sample field points in a vertical plane. We use
% this here to visualise the adiabatic energy levels + gravity in a slice
% plane through the potential.

RF = 3; % MHz
amp = 0.5 / 0.7; % Gauss
ap = AP.Calculator().CircularPolarised(RF, amp).DontUseParallel();
ap.Atom.F = 1;
disp(ap);

sampler = AP.Sampler.QuadrupolePlaneSampler(ap);
sampler.StartB = 2:0.2:4;
sampler.MeshIterations = 6;
sampler.ThetaRange = [ 0 pi*0.95 ];
sampler.RayNumber = 20;
sampler.Verbose = 1; sampler.QuadGrad = 100;
sampler.Sample();

E = sampler.E + Util.gpe(sampler.Z, 87);
tris = delaunay(sampler.X, sampler.Z);
t = trisurf(tris, sampler.X, sampler.Z, zeros(size(sampler.Y)), E(3,:));
xlabel('X (\mum)'); ylabel('Z (zmum)'); set(gcf, 'Color', 'w');
view([0,90]);
shading interp;

% We emphasise the field points and meshing by drawing faint edges on the
% triangles and black dots at the location of field points.
set(t, 'EdgeColor', 'k', 'EdgeAlpha', 0.1)
hold on; plot3(sampler.X, sampler.Z, ones(size(sampler.X)), '.k'); hold off;
