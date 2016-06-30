%% mrfSeries
% Define the rf amplitude sequences used in experiment and generate theory
% plots

rabi100 = [0.391 0.459 0.410]';

%% Series 2
% rabi100 = [0.391 0.459 0.410]';

c1 = 0.5*rabi100(1);
c2 = [0.20 0.40 0.60 0.80 0.85 0.90 0.95 1.00 1.05 1.10 1.20 1.30 1.40]*rabi100(2); % the final barrier heights
c3 = 1*rabi100(3);

s2minima = MRF.mrfPositions(c1, c2, c3);

figure(1)
plot(s2minima(:,1)/rabi100(2), s2minima(:,4),'.')

%% Series 3

c1 = 0.5*rabi100(1);
c2 = [1.4 1.3 1.2 1.1 1.0 0.9 0.8 0.5 0.4 0.3 0.2]*rabi100(2);... [1.4 1.3 1.2 1.1 1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2]*rabi100(2);
% had to omit 0.7 and 0.6 since didn't have a minimum here! Somthing must
% be wrong with this...
c3 = 1.11*rabi100(3);

s3minima = MRF.mrfPositions(c1, c2, c3);

figure(2)
plot(s3minima(:,1)/rabi100(2), s3minima(:,4),'.')

%% Series 1

c1 = 0.5*rabi100(1);
c2 = [0.4 0.6 0.8 1.0 0.9 1.1 1.2 1.3 1.4 1.5 1.6 0.5 0.7 0.3 0.85 0.95 1.05 1.15 0.2]*rabi100(2); % the final barrier heights
c3 = 1*rabi100(3);

s1minima = MRF.mrfPositions(c1, c2, c3);

figure(2)
plot(s1minima(:,1)/rabi100(2), s1minima(:,4),'.')

%% Testing

RFs = [3.0 3.6 4.2 ]';
c2 = [0.2]
% rabi100 = [0.391 0.459 0.410]'; % '100%' reference values measured in single-rf traps (MHz)

% Work with series 2 for now
% amp30 = c1 * rabi100(1);
% amp36 = c2 * rabi100(2); % the final barrier heights
% amp42 = c3 * rabi100(3);

qdrpGrad = 62.4511; % Gauss/cm at 20 A

Bs=(2.5:0.01:5); % range of Bs to consider 
ladno = 3;

rabis = [c1 c2(1) c3]';
[ F, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, rabis, 'iterations', 4, 'qdrpGrad', qdrpGrad);
F = MRF.sortEnergies(B, MRF.ladder(RFs, ladno, F));

figure(5)
plot(B,F)

lad = MRF.ladder(RFs, ladno, F);

figure(6)
plot(B, lad,'.')

figure(1)
plot(B, lad + repmat(MRF.gpe(B, qdrpGrad), size(lad,1), 1), '.');

Ftest = F(1,:);

Fworking = Ftest+repmat(MRF.gpe(B, qdrpGrad), size(Ftest,1), 1);

FworkingI = interp1(B, Fworking, Bs); % so this gives me equal spacing to plot lots of them!

figure(2)
plot(Bs, FworkingI,'.')
