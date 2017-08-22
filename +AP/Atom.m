classdef Atom < handle
    %ATOM Represents a species of atom.
    
    properties
        
        %F Hyperfine-spin of the atom.
        F = 1;
        
        %gFuB Magnetic dipole moment of the atom.
        gFuB = -0.7;
        
    end
    
    methods
       
        function a = Atom(F, gFuB)
           a.F = F;
           a.gFuB = gFuB;
        end
        
    end
    
end

