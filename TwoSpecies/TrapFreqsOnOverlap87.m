%% Trap frequencies on overlap
% We now attempt to determine the radial and axial trap frequencies for the
% MRF traps when the two clouds are overlapped.

%%
% Define trap parameters

QdrpGrad = 60;
zeemanSplit   = 2.9:0.2:3.4;
zeemanSplit   = 2.5:0.2:3.8;
gFuBB2_87 = 0.250; % 250 kHz \Omega_1
gFuBB1_87 = 0.4;   % 400 kHz \Omega_1
RFs = [ 4.5 3 ]';

%%
% Calculate the potential for 87 at the trapping location.

thetas=0:0.01:0.35;

xs = cell(1, length(thetas));
zs = cell(1, length(thetas));
vals = cell(1, length(thetas));

for i=1:length(thetas)
    theta = thetas(i);
    
     [ F, B ] = MRF.MeshedQuasiEnergies(zeemanSplit, RFs, [gFuBB1_87 gFuBB2_87]', 'iterations', 5, 'qdrpGrad', QdrpGrad, 'F', 1, 'theta', theta);
     F2 = MRF.sortEnergies(B, MRF.ladder(RFs, 30, F));
    
    % select a level; add these values to array
    trapped = 2;
    vals{i} = F2(trapped, :);
    
    plot(B, F2(trapped,:), '.'); pause(0.1);
    
    % calculate x and z:
    % convert B from MHz to Gauss
    B = B / 0.7; % gF uB in MHz/G
    x = ((1-cos(theta)) * (B / QdrpGrad).^2).^0.5 * 1e4;
    z = (cos(theta) * (B / QdrpGrad).^2 / 4).^0.5 * 1e4;
    
    xs{i} = x;
    zs{i} = z;
    
    fprintf('Calculating theta=%2.3f...\n', theta)
end

fprintf('Complete!\n')

%%
x = [];
z = [];
v = [];
for i=1:length(xs)
    x = [x xs{i}];
    z = [z zs{i}];
    v = [v vals{i}];
end

% Add in the gpe of points
v2 = v + gpe(-z, 87);

tri = delaunay(x,-z);
h = trisurf(tri, x, -z, v2);
hold on; plot3(x,-z,v2, 'k.'); hold off;
shading interp
view([0 90]);
xlabel('x (\mu m)');
ylabel('z (\mu m)');

box on

%%
% For each MRF 'slice', find the minimum position and it's energy.
minx = zeros(1, length(xs));
minz = zeros(1, length(xs));
minv = zeros(1, length(xs));
for i=1:length(xs)
   % within each slice
   energy = vals{i} + gpe(-zs{i}, 87);
   [~,j] = min(energy);
   temp = xs{i}; minx(i) = temp(j);
   temp = zs{i}; minz(i) = temp(j);
   minv(i) = energy(j);
end

hold on;
plot3(minx, -minz, minv, 'g.-');
hold off;

%%
% Walk along the bottom of this path to determine the 'arc length' s. This
% will be used for the harmoninc fit to determine the trap frequency.

s = zeros(1, length(minx));
for i=2:length(minx)
    s(i) = s(i-1) + ((minx(i)-minx(i-1)).^2 + (minz(i)-minz(i-1)).^2).^0.5;
end

plot(s, minv-minv(1));

freqs = Util.getTrapFreq(s, minv-minv(1), 87);
fprintf('The radial trap frequency is %.1f Hz.\n', freqs(2))

%% Time averaging
% Looks like the previous method does not give us nice parameters. Radial
% trap freqs are low, so the chemical potential too will be. Also, the 85
% confinement is weaker than that of the 87, which means it can just flow
% around the condensate rather than being forced to tunnel through it.

% radius of TOP orbit in um
TOPr = 150;

% We are only interested in relative energies. Trim off the large 'dc'
% component of potential that will otherwise lead to numerical imprecision.
v = v2 - min(v2(:));

% define interpolant potential
U = @(rq, zq) griddata(x, z, v, rq, zq);

% define line over which the TAAP will be sampled. First we must calculate
% the vertical minimum of the potential, which is shifted slightly due to
% the TAAP field.
zqs = minz(1) + (-15:1:15);
xqs = zeros(size(zqs));

% integrate points over 2*pi rotation of top field
Us = zeros(size(zqs));
phis = 0:2*pi/10:2*pi; phis(end) = [];
for phi=phis
    r = ((xqs - TOPr * cos(phi)).^2 + TOPr.^2 * sin(phi).^2).^0.5;
    
    Us = Us + U(r, zqs);
    
end
Us = Us / length(phis);

subplot(2,1,1);
plot(zqs, Us, '-'); hold on;

% Identify the new minimum energy location of the trap in the z-axis direction.
[~,i] = min(Us);
zqs   = zqs(i);

plot(zqs, Us(i), '.r'); hold off;
hold on; plot([1 1]*minz(1), ylim, ':r'); hold off;
title('Identification of vertical minimum');
xlabel('z (\mu m)');
ylabel('TAAP Potential (MHz)');

% determine sample line for radial direction trap freq.
xqs = 0:1:20;
zqs = zqs * ones(size(xqs));

% integrate points over 2*pi rotation of top field
Us = zeros(size(xqs));
phis = 0:2*pi/10:2*pi; phis(end) = [];
for phi=phis
    r = ((xqs - TOPr * cos(phi)).^2 + TOPr.^2 * sin(phi).^2).^0.5;
    U(r, zqs)
    Us = Us + U(r, zqs);
end
Us = Us / length(phis);

subplot(2,1,2);
plot(xqs(~isnan(Us)), Us(~isnan(Us)), '.-');% hold on; plot(xqs(isnan(vals)), 0, 'rx'); hold off;
freqs = Util.getTrapFreq(xqs(~isnan(Us)), Us(~isnan(Us)), 87);
title('Radial TAAP potential');
ylabel('TAAP Potential (MHz)');
xlabel('x (\mu m)');
fprintf('Trap frequency in TAAP = %.2f.\n', freqs(2))