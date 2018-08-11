classdef (Abstract) AbstractLineSampler < AP.Sampler.QuadrupoleSampler
    %ABSTRACTLINESAMPLER Samples the adiabatic potential along a line.
    %   Calculates eigenenergies and vectors of the dressed potential along
    %   a line between InitialSamplePoints(1) and (end). Employs meshing to
    %   refine the calculation, up to a specified number of iterations.
    %
    %   A quadrupole gradient can be specified, in which case the meshing
    %   will also account for gravitational sag.
    
    properties (SetAccess=protected, GetAccess=protected)
        
        %MB List of field points specified by position along line.
        mLambda = [];
        
        %ME List of field point eigenenergies, MHz, size(E,1)==atom.F
        mE = [];
        
        %MEV Eigenvectors. Rank 1,2 are eigenvectors, size(EV,3)==length(B)
        mEV = [];
        
    end
    
    properties (SetAccess=protected, GetAccess=public)
        
        %UNSORTEDEIGENENERGIES Get the unsorted eigenenergies
        UnsortedEigenenergies;
        
        %INITIALLAMBDA Initial lambda (positions) along line.
        InitialLambda = [];
        
    end
    
    properties (Dependent)
        
        %EIGENVECTORS Rank 1,2 are eigenvectors, size(EV,3) == length(B)
        Eigenvectors
        
        %ENERGIES Get the energies of this z-axis sampler at field points.
        Eigenenergies;
        
        Lambda
        
        %POTENTIALENERGIES Eigenenergy including gravitational sag
        PotentialEnergies;
       
    end
    
    properties
       
        %MESHITERATIONS Number of iterations to mesh calculation.
        MeshIterations = 3;
        
        %SORT Sort the eigenenergies by eigenvector similarity.
        Sort = 1;
        
        %VERBOSE Should the line sampler print debug information?
        Verbose = 1;
        
    end
    
    methods
        
        function instance = AbstractLineSampler(calculator)
            %LINESAMPLER3D Creates a new LineSampler3D instance.
            %   The constructor must be given a reference to an APCalculator
            %   object which will be used to sample the field points.
            instance = instance@AP.Sampler.QuadrupoleSampler(calculator);
        end
        
        Sample(instance);
        
        function coords = GetCoords(instance)
            %GETCOORDS Returns the coordinates at each field point.
            [B,theta,gamma] = lambdaToBTG(instance, instance.mLambda);
            coords = instance.TransformThetaGamma2XYZ(B, theta, gamma);
        end
        
        function Es = E(instance)
            %E Get the energies of this z-axis sampler at field points.
            
            if (instance.Dirty)
                instance.dirtyError()
            end
            
            Es = instance.mE;
        end
        
        function Esag = get.PotentialEnergies(instance)
           %POTENTIALENERGIES Get eigenenergies plus gravitational sag
            Esag = instance.Eigenenergies + Util.gpe(instance.Z, instance.APCalculator.Atom.Mass);
        end
        
        function Es = get.Eigenenergies(instance)
            %EIGENENERGIES Get the energies of this z-axis sampler at field points.
            Es = instance.E;
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
        
        function set.MeshIterations(instance, val)
            instance.MeshIterations = val;
            instance.Dirty = 1;
        end
        
        function set.Sort(instance, val)
            instance.Sort = val;
            instance.Dirty = 1;
        end
        
        function l = get.Lambda(instance)
            l = instance.mLambda;
        end
        
    end
    
    methods (Abstract)
       
        %lAMBDA2BTG Convert lambda to B, Theta, Gamma coordinates.
        %  B is Zeeman splitting in MHz
        [B,theta,gamma] = lambdaToBTG(instance, lambda);
        
        %ISHORIZONTAL Is the line horizontal? 
        h = IsHorizontal(instance)
        
    end
    
end