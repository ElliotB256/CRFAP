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
y = ap.GetDressedEnergies(x, pi/2, pi/2);
plot(x, y);