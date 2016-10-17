%%
% Attempt to render a slice through a shell in 3d. At first focus on a
% single RF, and sample the lower hemisphere.

RF = 3;
Rabi = 0.530;
qdrpGrad = 150;
Bs=(2.5:0.2:3.5);

thetas=0:0.05:(pi/2);
thetas = [thetas pi/2];

xs = cell(1, length(thetas));
zs = cell(1, length(thetas));
vals = cell(1, length(thetas));

for i=1:length(thetas)
    theta = thetas(i);
    
     [ F, B ] = MRF.MeshedQuasiEnergies(Bs, RF, Rabi, 'iterations', 3, 'qdrpGrad', qdrpGrad, 'F', 1, 'theta', theta);
     F2 = MRF.sortEnergies(B, MRF.ladder(RFs, 30, F));
    
    % select a level; add these values to array
    trapped = 2;
    vals{i} = F2(trapped, :);
    
    plot(B, F2(trapped,:), '.'); pause(0.1);
    
    % calculate x and z:
    % convert B from MHz to Gauss
    B = B / 0.7; % gF uB in MHz/G
    x = ((1-cos(theta)) * (B / qdrpGrad).^2).^0.5 * 1e4;
    z = (cos(theta) * (B / qdrpGrad).^2 / 4).^0.5 * 1e4;
%     if ~isreal(z)
%         z = -imag(z)
%     end
    % final factor = cm to microns
    
    xs{i} = x;
    zs{i} = z;
    
    fprintf('Calculating theta=%2.1f...\n', theta)
    x
end

fprintf('Complete!\n')


%%
% Triangulate the point cloud and present the potential surface

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
view([0 90]); axis equal
xlabel('x (\mu m)');
ylabel('z (\mu m)');
title(sprintf('$B''=%.0f$ G/cm, $\\Omega = %.2f$ MHz, $\\omega = %.2f$ MHz', qdrpGrad, Rabi, RF), 'Interpreter', 'Latex');

box on