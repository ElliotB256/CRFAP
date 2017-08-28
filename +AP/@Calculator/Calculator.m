classdef Calculator < handle
    %CALCULATOR Calculates dressed eigenstates and energies.
    %  The AP.Calculator class is used to calculate the eigenenergies and
    %  eigenstates of RF-dressed atoms in static magnetic fields. It is
    %  assumed that the underlying static field is a spherical quadrupole
    %  field.
    %
    %  The dressing rf field is assumed to be of the form:
    %   Brf = [ Bx sin(2pi RF t), By sin(2pi RF t + py), Bz sin(2pi RF t + pz)]
    %  where all field quantities are column vectors for multiple dressing
    %  fields. The quantisation axis is rotated from \hat{z} with respect
    %  to this field by an angle theta around the y axis, and then an angle
    %  gamma around the x' axis.
    
    properties (SetAccess=private)
        
        %RF Vector of dressing RFs in MHz
        RF;
        
        %BX Vector of field amplitudes in the x-direction (Gauss)
        BX;
        
        %BY Vector of field amplitudes in the y-direction (Gauss)
        BY;
        
        %BZ Vector of field amplitudes in the z-direction (Gauss)
        BZ;
        
        %PY Phase of Y-field (radians)
        PY;
        
        %PZ Phase of Z-field (radians)
        PZ;
        
        %PHASE Phase between different MRF components (radians)
        Phase;
        
        %ATOM Atomic species used for dressed energy calculation
        Atom;
        
        %PARALLEL Should the calculation occur using parallel loops?
        Parallel = 1;
        
    end
    
    methods
        
        function c = Calculator()
           c.OfSpecies('F', 1, 'OfSpecies', 87);
        end
        
        function context = CircularPolarised(context, RF, B, phase)
            %CIRCULARPOLARISED Configure for circular polarised dressing RF.
            %  Configures the dressing RF to be circularly polarised in the
            %  lab frame. The dressing field corresponds to a magnetic field
            %  of amplitude B (Gauss) that rotates around the z-axis at a
            %  frequency RF (MHz). These quantities can be vector quantities
            %  for MRF dressing.
            %
            %  The optional argument 'phase' allows the relative phase
            %  between dressing components to be set.
            %
            %  Syntax:
            %   ap = ap.CircularPolarised(RF, B)
            %   ap = ap.CircularPolarised(RF, B, phase)
            
            % Verify inputs
            if (nargin < 4)
                phase = zeros(size(RF));
            end
            
            % make columns
            RF = RF(:); B = B(:); phase = phase(:);
            
            if length(RF) ~= length(B) || length(RF) ~= length(phase)
                error('Vectors RF, B and phase must be of same length.');
            end
            
            % configure properties of the calculator.
            context.BZ = 0;
            context.PZ = 0;
            context.PY = pi/2;
            context.BX = B;
            context.BY = B;
            context.RF = RF;
            context.Phase = phase;
            
        end
        
        function result = IsCircularPolarised(context)
            %ISCIRCULARPOLARISED Is the system using circular polarised RF?
            
            % Requirement for circular polarised is that:
            %  'all phaseY = pi/2
            %  'all BX = BY
            %  'all Bz = 0
            
            result = all(context.PY == pi/2) && all(context.BX == context.BY) && all(context.BZ == 0);
            
        end
        
        function context = LinearPolarised(context, RF, B, phase)
            %LINEARPOLARISED Configure for linear polarised dressing RF.
            %  Configures the dressing RF to be linear polarised in the lab
            %  frame. The dressing field corresponds to a magnetic field of
            %  amplitude B (Gauss) that oscillates at a frequency RF (MHz)
            %  along the X axis. These quantities can be vector quantities
            %  for MRF dressing.
            %
            %  The optional argument 'phase' allows the relative phase
            %  between dressing components to be set.
            %
            %  Syntax:
            %   ap = ap.LinearPolarised(RF, B)
            %   ap = ap.LinearPolarised(RF, B, phase)
            
            % Verify inputs
            if (nargin < 4)
                phase = zeros(size(RF));
            end
            
            % make columns
            RF = RF(:); B = B(:); phase = phase(:);
            
            if length(RF) ~= length(B) || length(RF) ~= length(phase)
                error('Vectors RF, B and phase must be of same length.');
            end
            
            % configure properties of the calculator.
            context.BZ = 0;
            context.PZ = 0;
            context.PY = 0;
            context.BY = 0;
            context.BX = B;
            context.RF = RF;
            context.Phase = phase;
            
        end
        
        function result = IsLinearPolarised(context)
            %ISLINEARPOLARISED Is the system using linear polarised RF?
            
            % Requirement for linear polarised is that:
            %  'all BY = 0
            %  'all BZ = 0
            
            result = all(context.BY == 0) && all(context.BZ == 0);
            
        end
        
        function context = OfSpecies(context, varargin)
            %OFSPECIES Configure properties of the atomic Hamiltonian.
            
            ip = inputParser();
            ip.addParameter('OfSpecies', 87);
            ip.addParameter('F', 1);
            ip.parse(varargin{:});
            
            switch ip.Results.OfSpecies
                case 87
                    switch ip.Results.F
                        case 1
                            context.Atom = AP.Atom(1, -0.7);
                            return;
                    end
            end
            
            error('Unsupported atom species.');
        end
        
        function context = DontUseParallel(context)
            %DONTUSEPARALLEL Don't compute in parallel.
           context.Parallel = 0; 
        end
        
        function s = HilbertSpaceSize(context)
           %HILBERTSPACESIZE Gets the size of the Hilbert space.
            switch context.Atom.F
                case 1
                    s = 3;
                case 2
                    s = 5;
                otherwise
                    error('Unknown size of Hilbert Space');
            end
            
        end
        
    end
    
    methods
       
        %GETHAMILTONIAN Gets the rf-dressed Hamiltonian for calculations.
        H = GetHamiltonian(context, theta, gamma, omega0);
        
    end
    
end

