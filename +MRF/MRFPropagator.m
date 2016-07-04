function [ U ] = MRFPropagator( B, RF, BRF, time )
%PROPAGATOR Calculates the popagator operator that evolves a wavefunction
%for the specified time. This method uses Matlab's internal integrator
%ode45 which has a variable timestep. This method is faster than the more
%general 'Propagator', although it uses a more specialised Hamiltonian (for
%MRF).
% Syntax: MRFPropagator(B, RF, BRF, time)
%  Bs: energy splitting of the undressed Zeeman states in MHz.
%  RF: vector of dressing RFs (MHz)
%  BRF: vector of dressing RF Rabi frequencies (MHz)
%  time: time the propagator evolves the system for

    
    % We use that U*I = U, so evolve identity matrix for time to determine
    % propagator for that time.
    U = eye(3);
    [~,U] = ode45(@(t,ds) MRFPropTimeEv(t, ds, RF, BRF, B), [0 time], U(:));
    
    %recast U to sq matrix
    U = U(end,:);
    U = U([ 1 2 3; 4 5 6; 7 8 9;]');
end

