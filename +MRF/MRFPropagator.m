function [ U ] = MRFPropagator( B, RF, BRF, time )
%PROPAGATOR Calculates the MRF propagator for the given time. This method
%uses Matlab's internal integrator ode45 which has a variable timestep.
    
    % We use that U*I = U, so evolve identity matrix for time to determine
    % propagator for that time.
    U = eye(3);
    [~,U] = ode45(@(t,ds) MRFPropTimeEv(t, ds, RF, BRF, B), [0 time], U(:));
    
    %recast U to sq matrix
    U = U(end,:);
    U = U([ 1 2 3; 4 5 6; 7 8 9;]');
end

