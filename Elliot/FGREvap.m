%% Evap via FGR

%% Single RF
% We represent the atom/photon states in the basis of |m> (x) |n_1>. We
% truncate the maximum number of photons in mode 1 to N.
%
% Basis looks like this:
%  | mF, n_1 >
%  | -1, 0 >
%  | 0 , 0 >
%  | 1 , 0 >
%  | -1, 1 >
%  | 0 , 1 >
%  | 1 , 1 >
%  ....



N = 20;
rf = 1;%2*pi*1;

H_a = diag( [ -1 0 1 ]');
%H_a = diag(repmat([ -1 0 1 ]', N, 1));
H_rf1 = diag( rf * (1:N));

% identity matrices:
Ia      = eye(3);
Irf1    = eye(N);

H_nonInt = kron(Irf1, H_a) + kron(H_rf1, Ia);

% Add in interaction terms:
g = 0.5;
% use off-diagonal indices to construct interaction matrix
H_int   = zeros(size(H_nonInt,1));
i = 1:3:size(H_nonInt,1)-4;
H_int(sub2ind(size(H_int),i, i+4))   = g;
H_int(sub2ind(size(H_int),i+1, i+5)) = g;
H_int(sub2ind(size(H_int),i+4, i))   = g;
H_int(sub2ind(size(H_int),i+5, i+1)) = g;
imagesc(H_int);

label = labelStates( [ 3 N ] );

% We now can formulate the complete Hamiltonian
H = H_nonInt + H_int;
imagesc(H);

% Calculate the eigenstates of this Hamiltonian
[E,~] = eig(H);
imagesc(E); title('Eigenvalue decomposition');

% Define transformation from 'bare' basis to dressed state basis
toDressedBasis = @(U) E \ U * E;


%% Adding in the probe RF
% Now for the fun part: We need to mix in a weak probe RF. We expand the
% basis of the Hamiltonian to now include photo energies of the weak probe
% RF. For sake of computational power we limit the number of photons in the
% second probe RF to a lower number - this is fine for a weak probe.

prbRF = 1.5;
gPrb = 0.01;

Nprb = 3;
IrfPrb = eye(Nprb);
% self-energy of probe photons
H_rfPrb = diag(prbRF * (1:Nprb));

% Generate the probe RF matrix. The interaction will lower both m and n2 or
% raise both m and n2. I'll make a quick and dirty enumeration for now to
% get the general form.
H_int_prb = zeros(prod([ 3 N Nprb ]));

% occupation numbers for the bare states are easy to use.
labels = labelStates( [ 3 N Nprb ] );

for i=1:size(H_int_prb, 1)
    for j=1:size(H_int_prb, 2)
        
        % get labels for the current two states.
        a = labels(i,:);
        b = labels(j,:);
        
        % set H_int_rf to non-zero where we are raising or lowering both
        % atomic level and rf2 level.
        if a(3) == b(3) + 1 && a(1) == b(1) + 1
            H_int_prb(i,j) = 1;
        end
        
        if a(3) == b(3) - 1 && a(1) == b(1) - 1
            H_int_prb(i,j) = 1;
        end
        
    end
end

% multiply by probe RF strength
H_int_prb = H_int_prb * gPrb;

%% Calculating FGR rates
% We now have a representation of the probe RF interaction in the basis
% |mF,N1,prbN>. Let's get the eigenstates of the RF1 dressed atoms and
% express them in this basis too.

vec = zeros(size(E,1),1);
DS = kron(IrfPrb, E);

% Transform the probe rf interaction into this dressed state representation
probeIntDS = DS \ H_int_prb * DS;

% Work out matrix elements for transitions from one chosen dressed state to
% all others. Enumerate through final states and fill 'rates' with matrix
% element between 'vec' and each final state.

% initial state is from a selected eigenstate of the single RF dressed
% hamiltonian (trapped state) and with nPrb = 1 so we can have both probe
% absorption and emission events.
initial = vec; initial(1) = 1; initial = kron([ 0 1 0 ]', E * initial);

rates = [];
% final states are all other states BUT NOT THIS ONE!
% i indexes over each dressed state. j indexes over N2 states
for j=[1 2 3];
    for i=1:size(E,1)
        final = vec; final(j) = 1; final = kron([ 0 1 0 ]', E * final);
        rates(end+1) = final' * probeIntDS * initial;
    end
end