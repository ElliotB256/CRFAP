%%
RFs = [3 3.6 4.2 ]';
BRFs = 0.4 * [ 0.5 0.2 1.1 ]';
eigF = [];
periodicity = 10*pi/3;
Bs=(2.5:0.05:5.5);

for i=1:length(Bs)
    B = Bs(i);
    %H = MRF.Hamiltonian(B, RFs, BRFs);
    %U = MRF.propagator(H, periodicity);
    U = MRF.MRFPropagator(B, RFs, BRFs, periodicity);
    eigF(:,i) = sort(angle(eig(U)))/periodicity;
end

plot(Bs, eigF','xr'); hold on;

%%
% With more general propagator notation
RFs = [3 3.6 4.2 ]';
BRFs = 0.8 * [ 0.5 0.2 1.1 ]';
eigF2 = [];
periodicity = 10*pi/3;
deltaB = 0.05;
Bs=(2.5:deltaB:5.5);

% Note: At present, the Bs are specified in MHz. This means the factor gF
% uB has been incoporated into the B specified here. Think of it more as
% the energy splitting of the bare states.

for i=1:length(Bs)
    B = Bs(i);
    H = MRF.Hamiltonian(B, RFs, BRFs);
    U = MRF.Propagator(H, periodicity);
    eigF2(:,i) = sort(angle(eig(U)))/periodicity;
end

hold on; plot(Bs, eigF2', '.'); hold off;

% repmat and produce level structure...
gcd = 2*pi/periodicity;
offsets = gcd*(-3:1:3)';
lad = repmat(eigF2, length(offsets), 1) + repmat(offsets, 3, size(eigF2, 2)); 
plot(Bs,lad')

% mesh in finer detail around stationary points. First, find stationary
% points
d = diff(eigF2, 1, 2)';
d2 = diff(sign(d), 1, 1);
i = find(abs(d2(:,1)) > 0.5);

% Having identified these stationary points, mesh them in finer detail.
cB = Bs(i); %centres to mesh
MeshN = 10;
mesh = -deltaB:2*deltaB/MeshN:deltaB;
% remove duplicate points from mesh
mesh = mesh(2:end-1);
mesh = mesh([1:(MeshN/2 -1) (MeshN/2+1):size(mesh,2)]);

% define new field points
Bs2 = repmat(cB, size(mesh, 2), 1) + repmat(mesh', 1, size(cB, 2));
Bs2 = Bs2(:)';

eigF3 = [];

% Calculate energies
for i=1:length(Bs2)
    B = Bs2(i);
    H = MRF.Hamiltonian(B, RFs, BRFs);
    U = MRF.Propagator(H, periodicity);
    eigF3(:,i) = sort(angle(eig(U)))/periodicity;
end

% plot old data and new data remeshed over the top
B = [Bs Bs2];
E = [eigF2 eigF3];
