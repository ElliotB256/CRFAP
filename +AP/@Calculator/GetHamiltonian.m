function [ H ] = GetHamiltonian( context, omega0, theta, gamma )
%GETHAMILTONIAN Gets the rf-dressed Hamiltonian for calculations.
%  Creates a function handle H(t) that returns the Hamiltonian of the
%  rf-dressed system at a time t.
% 
%  The Hamiltonian is represented as a matrix in terms of the basis states
%  |m_F>. The order is such that |state>(0) = amplitude of most negative mF
%  state. For example, for F=1 the vector [1,0,0] represents the state
%  |m_F=-1>, the vector [0,1,0] represents |m_F=0> and the vector [0,0,1]
%  represents |m_F=1>.
%  
%  The rotation of the quantisation axis by angles theta, gamma is
%  described in the description for AP.Calculator.
%  
%  Note that omega0 is equal to gFuBB0, where B0 is the static field. For
%  atoms with sign(gF)<0 this quantity is negative. For this function,
%  omega0 should always be specified as positive, and sign(gF)abs(omega0)
%  will be passed to the Hamiltonian sub-functions.
%  
%  Automatically attempts to apply approximations to simplify the terms
%  required for the Hamiltonian.
% 
%  Syntax: ap.GetHamiltonian(omega0, theta, gamma)

% special case for atoms along the z-axis.
bottomOfShell = theta == 0;

% Unpack class properties into local variables. This will improve execution
% speed of the anonymous functions.
BX = context.BX(:); 
BY = context.BY(:);
BZ = context.BZ(:);
PY = context.PY(:);
PZ = context.PZ(:);
RF = context.RF(:);
phase = context.Phase(:);
gFuB = context.Atom.gFuB;
omega0 = abs(omega0)*sign(gFuB);
F = context.Atom.F;

switch F
   
    case 1
        
        if context.IsCircularPolarised()
            gFuBB = gFuB * BX;
            if (bottomOfShell)
                H = @(t) AP.Hamiltonian.F1CircPolBottomOfShell(t, omega0, RF, gFuBB, phase); 
                return;
%             else
%                 H = @(t) AP.Hamiltonian.F1CircPol(t, omega0, RF, gFuBB, theta, phase); 
%                 return;
            end
            
        elseif context.IsLinearPolarised()
            gFuBB = gFuB * BX;
            if (bottomOfShell)
                H = @(t) AP.Hamiltonian.F1LinPolBottomOfShell(t, omega0, RF, gFuBB, phase );
                return;
            else
                H = @(t) AP.Hamiltonian.F1LinPol( t, omega0, RF, theta, gamma, gFuBB, phase );
                return;
            end
        end
        
        H = @(t) AP.Hamiltonian.F1General( t, omega0, RF, theta, gamma, gFuB*BX, gFuB*BY, gFuB*BZ, PY, PZ, phase );
    
    otherwise    
        if (bottomOfShell)
            H = @(t) AP.Hamiltonian.BottomOfShell( t, F, omega0, RF, gFuB*BX, gFuB*BY, gFuB*BZ, PY, PZ, phase );
        else
            H = @(t) AP.Hamiltonian.General( t, F, omega0, RF, theta, gamma, gFuB*BX, gFuB*BY, gFuB*BZ, PY, PZ, phase );
        end
end

end

