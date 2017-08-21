function [ H ] = GetHamiltonian( context, omega0, theta, gamma )
%GETHAMILTONIAN Gets the rf-dressed Hamiltonian for calculations.
%  Creates a function handle H(t) that returns the Hamiltonian of the
%  rf-dressed system at a time t.
%  
%  Automatically attempts to apply approximations to simplify the terms
%  required for the Hamiltonian.
% 
%  Syntax: ap.GetHamiltonian(omega0, theta, gamma)

% special case for atoms along the z-axis.
bottomOfShell = theta == 0;

% Unpack class properties into local variables. This will improve execution
% speed of the anonymous functions.
BX = context.BX; 
BY = context.BY;
BZ = context.BZ;
PY = context.PY;
PZ = context.PZ;
RF = context.RF;
phase = context.Phase;
gFuB = context.gFuB;

switch context.F
   
    case 1
        
        if context.IsCircularPolarised()
            gFuBB = gFuB * BX;
            if (bottomOfShell)
                H = @(t) F1CircPolBottomOfShell(t, omega0, RF, gFuBB, phase); 
                return;
            else
                H = @(t) F1CircPol(t, omega0, RF, gFuBB, theta, phase); 
                return;
            end
            
        elseif context.IsLinearPolarised()
            gFuBB = gFuB * BX;
            if (bottomOfShell)
                H = @(t) F1LinPolBottomOfShell(t, omega0, RF, gFuBB, phase );
                return;
            else
                H = @(t) F1LinPol( t, omega0, RF, theta, gamma, gFuBB, phase );
                return;
            end
        end
        
        error('Not implemented for this system.');
    case 2
    
        error('Not implemented for this system.');
end

end

