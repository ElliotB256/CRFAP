    function eigF = compute(B, RF, gFuBB, p, periodicity)
    %COMPUTE Computes MRF eigenenergies. Used in GetQuasiEnergies.
    % B: zeeman splitting, MHz
    % RF: dressing frequencies, MHz
    % gFuBB: the value of gF * Bohr magneton * field amplitude B. equal to
    %        rabi freq for circ polarisations.
    % p: input argument struct from GetQuasiEnergies
    % periodicity: periodicity of the MRF system.
    
        cH = MRF.Hamiltonian(B, RF, gFuBB, ...
            'F', p.Results.F, ...
            'theta', p.Results.theta, ...
            'phase', p.Results.phase, ...
            'polarisation', p.Results.polarisation);
        cU = MRF.Propagator(cH, periodicity, p.Results.F);
        eigF = sort(angle(eig(cU)))/periodicity;
    end