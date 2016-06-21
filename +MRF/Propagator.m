function [ U ] = Propagator( H, time )
%PROPAGATOR Calculates the propagator for the given Hamiltonian after
%specified time. This method takes longer than the MRF propagator but is
%more general, allowing any Hamiltonian to be specified.
    
    % We use that U*I = U, so evolve identity matrix for time to determine
    % propagator for that time.
    U = eye(3);
    [~,U] = ode45(@(t,ds) PropTimeEv(t, ds, H), [0 time], U(:));
    
    %recast U to sq matrix
    U = U(end,:);
    U = U([ 1 2 3; 4 5 6; 7 8 9;]');
end

