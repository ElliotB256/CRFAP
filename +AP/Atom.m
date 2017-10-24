classdef Atom < handle
    %ATOM Represents a species of atom.
    
    properties
        
        %F Hyperfine-spin of the atom.
        F = 1;
        
        %gFuB Magnetic dipole moment of the atom.
        gFuB = -0.7;
        
        %MASS Mass of atom, in atomic mass units
        Mass = 87;
        
    end
    
    methods
        
        function a = Atom(F, gFuB, mass)
            if nargin <3
                mass = 87;
            end
            a.F = F;
            a.gFuB = gFuB;
            a.Mass = mass;
        end
        
    end
    
end

