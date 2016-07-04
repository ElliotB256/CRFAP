function [ U ] = propagator2( B, RF, BRF, time, prec )
%PROPAGATOR Calculates the propagator for the given Hamiltonian for
%given time. H is a struct describing the Hamiltonian    
    
    % time integration
    dt = prec * time;
    
    % We use that U*I = U, so evolve identity matrix for time to determine
    % propagator for that time.
    U = eye(3);
    for t=0:dt:time
        h = MRF.calcH(B, RF, BRF, t);
        deltaM = -1i * dt * ( h * U );
        U = U + deltaM;
    end
end

