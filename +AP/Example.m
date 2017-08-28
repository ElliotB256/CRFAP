% Configure AP Calculator objects and calculate the dressed energy
% levels.

%% Circular polarised

RF = 3; % MHz
amp = 0.5; % Gauss
ap = AP.Calculator().CircularPolarised(RF, amp);
disp(ap);

x = 2:0.1:4;
y = ap.GetDressedEnergies(x);
plot(x, y);

%% Linear polarised

RF = 3; % MHz
amp = 0.5; % Gauss
ap = AP.Calculator().LinearPolarised(RF, amp);
disp(ap);

x = 2:0.1:4;
y = ap.GetDressedEnergies(x);
plot(x, y);

%% Linear polarised - general rotation

RF = 3; % MHz
amp = 0.5; % Gauss
ap = AP.Calculator().LinearPolarised(RF, amp);
disp(ap);

x = 2:0.1:4;

%% Tests
% Check local quantisation axis parallel to x (and therefore to B):
y = ap.GetDressedEnergies(x, pi/2, 0);
plot(x, y);

%%
% Check local quantisation axis parallel to y (perpendicular to B):
y = ap.GetDressedEnergies(x, pi/2, pi/2);
plot(x, y);

%%
% Check local quantisation axis parallel to z (perpendicular to B):
y = ap.GetDressedEnergies(x, 0.001, 0);
plot(x, y);

%%
% Use of ZAxis sampler
RF = 3; % MHz
amp = 0.5 / 0.7; % Gauss
ap = AP.Calculator().CircularPolarised(RF, amp).DontUseParallel();
ap.Atom.F = 2;
disp(ap);

sampler = AP.Sampler.ZAxisSampler(ap);
sampler.StartB = 2:0.2:4; sampler.Sort = 1;
sampler.MeshIterations = 8;
sampler.Verbose = 1; sampler.QuadGrad = 100;
sampler.Sample();
plot(sampler.Z, sampler.E, '.-');

%%
% More complicated example - linear for F=1, F=2 of Rb87
RF = 3; % MHz
amp = 0.5 / 0.7; % Gauss
% 
ap = AP.Calculator().LinearPolarised(RF, amp).DontUseParallel().OfSpecies('species', 87, 'F', 1);
ap.PY = pi/2; ap.BY = ap.BX * 0.2;
sampler = AP.Sampler.ZAxisSampler(ap);
sampler.StartB = 2:0.2:4; sampler.Sort = 1;
sampler.MeshIterations = 8; sampler.Verbose = 1;
sampler.Sample();
clf; plot(sampler.B, sampler.E, '.-r'); hold on;

ap = AP.Calculator().LinearPolarised(RF, amp).DontUseParallel().OfSpecies('species', 87, 'F', 2);
ap.PY = pi/2; ap.BY = ap.BX * 0.2;
sampler = AP.Sampler.ZAxisSampler(ap);
sampler.StartB = 2:0.2:4; sampler.Sort = 1;
sampler.MeshIterations = 8; sampler.Verbose = 1;
sampler.Sample();
plot(sampler.B, sampler.E, '.-k'); hold off;
pause(1);