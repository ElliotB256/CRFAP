%% Evap via FGR

%% Single RF
% We represent the atom/photon states in the basis of |m> (x) |n_1>. We
% truncate the maximum number of photons in mode 1 to N.
%
% Basis looks like this:
%  | mF, n_1, n_3 >
%  | -1, 0  , 0 >
%  | 0 , 0  , 0 >
%  | 1 , 0  , 0 >
%  | -1, 1  , 0 >
%  | 0 , 1  , 0 >
%  | 1 , 1  , 0 >
%  ....
%
% where mF=atomic state, n_1 = dressing RF, n_2 = evap/probe


N = 10;
Np = 3;
RF = 1; % * 2 * pi
RFp = 0.5;

gRF = 0.5;
gRFp = 0.01;

H_a = diag( [ -1 0 1 ]');
H_RF = diag( RF * (1:N));
H_RFp = diag( RFp * (1:Np));

% identity matrices:
I_a      = eye(3);
I_RF     = eye(N);
I_RFp    = eye(Np);

% define basis: occupation numbers for the bare states are easy to use.
states = labelStates( [ 3 N Np ] );
merge = @(a,b,c) kron(I_RFp, kron(I_RF, a)) + kron(I_RFp, kron(b, I_a)) + kron(c, kron(I_RF, I_a));

% Calculate the Hamiltonian of the atomic energy plus the field energies.
H0 = merge(H_a, H_RF, H_RFp);

H_int = zeros(size(states, 1));
% Construct interaction matrix in a more general way using labelled states
for i=1:size(H_int,1)
    for j=1:size(H_int,2)
        
        a = states(i, :);
        b = states(j, :);
        
        % both are raised with no effect on probe N:
        if a(1) == b(1) + 1 && a(2) == b(2) + 1 && a(3) == b(3)
            H_int(i,j) = 1;
        end
        
        % both are lowered with no effect on probe N:
        if a(1) == b(1) - 1 && a(2) == b(2) - 1 && a(3) == b(3)
            H_int(i,j) = 1;
        end
    end
end
H_int = H_int * gRF;

% We now can formulate the complete Hamiltonian
H = H0 + H_int;
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

% Generate the matrix for the interaction of the probe RF. The interaction
% will lower both m and n2 or raise both m and n2.
H_probeInt = zeros(size(states, 1));
for i=1:size(H_probeInt, 1)
    for j=1:size(H_probeInt, 2)
        
        % get labels for the current two states.
        a = states(i,:);
        b = states(j,:);
        
        % both mF and probe photon number are raised
        if a(3) == b(3) + 1 && a(1) == b(1) + 1
            H_probeInt(i,j) = 1;
        end
        
        % both mF and probe photon number are lowered
        if a(3) == b(3) - 1 && a(1) == b(1) - 1
            H_probeInt(i,j) = 1;
        end
        
    end
end
H_probeInt = H_probeInt * gRFp;

%% Calculating FGR rates
% We now have a representation of the probe RF interaction in the basis
% |mF,N1,prbN>. Let's transform it into the basis of the dressed states.

% Transform the probe rf interaction into this dressed state representation
probeIntDS = E \ H_probeInt * E;

% Work out matrix elements for transitions from one chosen dressed state to
% all others. Enumerate through final states and fill 'rates' with matrix
% element between 'vec' and each final state.

vec = zeros(size(E,1),1);

% initial state is from a selected eigenstate of the single RF dressed
% hamiltonian (trapped state) and with nPrb = 1 so we can have both probe
% absorption and emission events.
initial = vec; initial(1) = 1; initial = E * initial;

rates = [];
% final states are all other states BUT NOT THIS ONE!
% i indexes over each dressed state. j indexes over N2 states
for j=[1 2 3];
    for i=1:size(E,1)
        final = vec; final(i) = 1; final = E * final;
        rates(end+1) = conj(final)' * probeIntDS * initial;
        %rates(end+1) = conj(final)' * initial;
    end
end

totalRate = sum(abs(rates(:)).^2);