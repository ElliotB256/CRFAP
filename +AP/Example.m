% Configure an AP Calculator object and calculate the dressed energy
% levels.

RF = 3; % MHz
amp = 0.5; % Gauss
ap = AP.Calculator().CircularPolarised(RF, amp);
disp(ap);

x = 2:0.1:4;
y = ap.GetDressedEnergies(x);
plot(x, y);