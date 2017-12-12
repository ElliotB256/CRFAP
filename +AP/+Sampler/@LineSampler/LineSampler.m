classdef LineSampler < AP.Sampler.QuadrupoleSampler
    %LINESAMPLER Samples the adiabatic potential along a line.
    %   Calculates eigenenergies and vectors of the dressed potential along
    %   an axis of constant \theta, \gamma. Employs meshing to refine the
    %   calculation, up to a specified number of iterations.
    %
    %   Field point coordinates are stored internally in terms of the
    %   zeeman splitting, but can be requested as coordinates.
    %
    %   A quadrupole gradient can be specified, in which case the meshing
    %   will also account for gravitational sag.
    
    properties (SetAccess=protected, GetAccess=protected)
        
        %MB List of field points zeeman splittings, MHz.
        mB = [];
        
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
        
        %B Get the field coordinates of this z-axis sampler (MHz).
        B;
        
        %POTENTIALENERGIES Eigenenergy including gravitational sag
        PotentialEnergies;
        
    end
    
    properties
        
        %STARTB List of field point coordinates, Zeeman splitting in MHz.
        StartB = [];
        
        %MESHITERATIONS Number of iterations to mesh calculation.
        MeshIterations = 3;
        
        %SORT Sort the eigenenergies by eigenvector similarity.
        Sort = 1;
        
        %VERBOSE Should the line sampler print debug information?
        Verbose = 1;
        
        %THETA Constant theta that describes this line.
        Theta;
        
        %GAMMA Constant gamma that describes this line.
        Gamma;
        
    end
    
    methods
        
        function instance = LineSampler(calculator, theta, gamma)
            %LINESAMPLER Creates a new LineSampler instance.
            %   The constructor must be given a reference to an APCalculator
            %   object which will be used to sample the field points.
            
            if nargin < 3
                gamma = 0;
            end
            
            if nargin < 2
                theta = pi;
            end
            
            instance = instance@AP.Sampler.QuadrupoleSampler(calculator);
            instance.Theta = theta;
            instance.Gamma = gamma;
        end
        
        Sample(instance);
        
        function coords = GetCoords(instance)
            %GETCOORDS Returns the coordinates at each field point.
            coords = instance.TransformThetaGamma2XYZ(instance.mB, instance.Theta, instance.Gamma);
        end
        
        function Bs = get.B(instance)
            %B Get the field coordinates of this z-axis sampler.
            %   Returns zeeman splitting in MHz.
            
            if (instance.Dirty)
                instance.dirtyError()
            end
            
            Bs = instance.mB;
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
        
        function Es = Eigenstates(instance)
            %EIGENSTATES Get eigenstates of the dressed rf calculation.
            
            if (instance.Dirty)
                instance.dirtyError()
            end
            
            Es = instance.mE;
        end
        
        function vectors = get.Eigenvectors(instance)
            %EIGENVECTORS Get eigenvectors from this linesampler.
            
            if (instance.Dirty)
                instance.dirtyError()
            end
            
            vectors = instance.mEV;
        end
        
        function set.StartB(instance, val)
            instance.StartB = val;
            instance.Dirty = 1;
        end
        
        function set.MeshIterations(instance, val)
            instance.MeshIterations = val;
            instance.Dirty = 1;
        end
        
        function set.Sort(instance, val)
            instance.Sort = val;
            instance.Dirty = 1;
        end
        
        function set.Theta(instance, val)
            instance.Theta = val;
            instance.Dirty = 1;
        end
        
        function set.Gamma(instance, val)
            instance.Gamma = val;
            instance.Dirty = 1;
        end
        
    end
    
end

