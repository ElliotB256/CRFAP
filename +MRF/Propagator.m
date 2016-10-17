function [ U ] = Propagator( H, time, F )
%PROPAGATOR Calculates the propagator operator that evolves a wavefunction
%for the specified time. This method uses Matlab's internal integrator
%ode45 which has a variable timestep.
% Syntax: MRFPropagator(H, time)
%  H: The Hamiltonian describing the system
%  time: time the propagator evolves the system for

if nargin < 3
    F = 1;
end

% We use that U*I = U, so evolve identity matrix for time to determine
% propagator for that time.
if F == 1
    U = eye(3);
    [~,U] = ode45(@(t,ds) PropTimeEv(t, ds, H), [0 time], U(:), odeset('RelTol', 1e-5));
    %recast U to sq matrix
    U = U(end,:);
    U = U([ 1 2 3; 4 5 6; 7 8 9;]');
elseif F == 2
    U = eye(5);
    [~,U] = ode45(@(t,ds) PropTimeEvF2(t, ds, H), [0 time], U(:), odeset('RelTol', 1e-3));
    %recast U to sq matrix
    U = U(end,:);
    U = U([1:5; 6:10; 11:15; 16:20; 21:25]');
else
    error('Invalid F: must be 1 or 2.');
end

end