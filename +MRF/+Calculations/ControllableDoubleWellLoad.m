%% Controllable double well loading
% Calculating theoretical estimates of the percentage of atoms loaded into
% upper well as a function of barrier height.
% Elliot

% First up: Let's look at the extreme limit of a fully submerged barrier.
% In this limit the barrier is lower than thermal energy, allowing
% thermalisation of the cloud between the two wells. In this limit the
% proportion of atoms between the two wells is a function of the initial
% temperature only and if truly adiabatic it should be independent of the
% final barrier height.

initialT = 1e-6;
K2MHz = @(x) x * Constants.kB / Constants.h;

% Construct the double well potential

Bs = 2.5:0.025:4.5;
BaseRabi = [ 338.7 401.2 357.0]' * 1e-3;
RFs = [3.0 3.6 4.2 ]';
qdrpGrad = 5*62.4511*0.96; % Gauss/cm at 100 A - check linearity cf 20A?


Rabi = BaseRabi .* [ 0.68 0.9 1.0 ]';

% revisited on 01/09/2016
Rabi = [ 0.2616 0.4512*0.8 0.4426 ]';
qdrpGrad = 59.5683;

[ F, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, Rabi, 'iterations', 0, 'qdrpGrad', qdrpGrad);
F = MRF.sortEnergies(B, MRF.ladder(RFs, 3, F));
lad = MRF.ladder(RFs, 3, F);
lad = lad + repmat(MRF.gpe(B, qdrpGrad), size(lad,1), 1);
plot(B, lad, '-');

%%
% pick out a trapped level

ti = 1;
trapped = lad(ti, :);

% Plot thermal probability distribution in this potential. Zero energy at
% minimum of the trapped level.

minE = min(trapped(:));
trapped = trapped - minE;

% Create partition function for ensemble
p = @(E) -E*1e6*Constants.h/(Constants.kB * initialT);

% Note: probably need to consider degeneracy due to available momentum
% states?
Z = sum(exp(1).^p(trapped(:)));
pDist = exp(1).^p(trapped(:)) / Z;

plot(B, trapped, '.-'); hold on;
plot(B, pDist); hold off; ylim([0 max(pDist)]);

% It is going to be very difficult to analyse the loading based on varying
% the barrier height as the cloud won't be thermalised over the barrier for
% all such heights.
% 
% Maybe a better approach would be to actually completely lower the
% barrier, then raise it again, for different RF1 amplitudes? This will
% raise/lower the height of the rf 1 well wrt the RF3 well. For high
% RF2 amplitudes it will be possible to analyse this just assuming
% thermalisation of the ensemble between the two wells.