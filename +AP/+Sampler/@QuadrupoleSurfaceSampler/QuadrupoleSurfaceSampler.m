classdef QuadrupoleSurfaceSampler < AP.Sampler.QuadrupoleSampler
    %QUADRUPOLESURFACESAMPLER Samplers points on resonant spheroid
    
    properties
        
        %THETA Vertical angles alpha at which to sample
        Theta;
        
        %GAMMA Radial angles at which to sample
        Gamma;
        
        %RF RF resonance to sample
        RF;
        
        %VERBOSE Should the line sampler print debug information?
        Verbose = 1;
        
    end
    
    properties (SetAccess=protected, GetAccess=protected)
        
        %ME List of field point eigenenergies, MHz, size(E,1)==atom.F
        mE = [];
        
        %MEV Eigenvectors. Rank 1,2 are eigenvectors, size(EV,3)==length(B)
        mEV = [];
        
    end
    
    properties (Dependent)
        
        %EIGENVECTORS Rank 1,2 are eigenvectors, size(EV,3) == length(B)
        Eigenvectors
        
        %ENERGIES Get the energies of this z-axis sampler at field points.
        Eigenenergies;
        
        %POTENTIALENERGIES Eigenenergy including gravitational sag
        PotentialEnergies;
        
    end
    
    methods
        
        function instance = QuadrupoleSurfaceSampler(calculator, theta, gamma)
            %QUADRUPOLESURFACESAMPLER Creates a new QuadrupoleSurfaceSampler instance.
            %   The constructor must be given a reference to an APCalculator
            %   object which will be used to sample the field points.
            
            if nargin < 3
                gamma = linspace(0, 2 * pi, 12);
                gamma = gamma(1:end-1);
            end
            
            if nargin < 2
                theta = linspace(0, pi/2, 6);
            end
            
            instance = instance@AP.Sampler.QuadrupoleSampler(calculator);
            instance.Theta = theta;
            instance.Gamma = gamma;
            instance.RF = calculator.RF(1);
        end
        
        function thetas = GetThetas(instance)
            a = instance.Theta(:);
            b = instance.Gamma(:);
            thetas = repmat(a', length(b), 1); thetas = thetas(:);
        end
        
        function gammas = GetGammas(instance)
            a = instance.Theta(:);
            b = instance.Gamma(:);
            gammas = repmat(b, 1, length(a)); gammas = gammas(:);
        end
        
        function Es = E(instance)
            %E Get the energies of this z-axis sampler at field points.
            
            if (instance.Dirty)
                instance.dirtyError()
            end
            
            Es = instance.mE;
        end
        
        function Es = get.Eigenenergies(instance)
            %EIGENENERGIES Get the energies of this z-axis sampler at field points.
            Es = instance.E;
        end
        
        function Esag = get.PotentialEnergies(instance)
           %POTENTIALENERGIES Get eigenenergies plus gravitational sag
            Esag = instance.Eigenenergies + Util.gpe(instance.Z, instance.APCalculator.Atom.Mass);
        end
        
        function coords = GetCoords(instance)
            %GETCOORDS Returns the coordinates at each field point.
            t = instance.GetThetas();
            g = instance.GetGammas();
            coords = instance.TransformThetaGamma2XYZ(instance.RF, t, g);
        end
        
        Sample(instance);
        
    end
    
end

