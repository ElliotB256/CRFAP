function [ U ] = EulerPropagator( H, time )
%PROPAGATOR Calculates the propagator for the given Hamiltonian for
%given time. This is general and works for any propagator, but uses a
%somewhat primitive Euler integration. The time integration was found to be
%significantly slowed down by the use of a function handle for the
%Hamiltonian, which added a factor of 10 slowdown.
    
    % time integration
    dt = 0.0001 * time;
    
    % We use that U*I = U, so evolve identity matrix for time to determine
    % propagator for that time.
    U = eye(3);
    for t=0:dt:time
        h = H(t);
        deltaM = -1i * dt * ( h * U );
        U = U + deltaM;
    end
end

