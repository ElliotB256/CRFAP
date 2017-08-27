classdef ZAxisSampler < AP.Sampler.AbstractSampler & matlab.mixin.SetGet
    %ZAXISSAMPLER Samples the adiabatic potential along the z axis.
    %   Calculates eigenenergies and vectors of the dressed potential along
    %   the vertical axis below the quadrupole. Employs meshing to refine
    %   the calculation, up to a specified number of iterations.
    %
    %   Field point coordinates are stored internally in terms of the
    %   zeeman splitting, but can be requested in Z coordinates.
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
        
        %APCALCULATOR Reference to APCalculator used for field point calc.
        APCalculator;
        
    end
    
    properties (Transient)
        %DIRTY The instance is dirty if properties are changed/before calc.
        Dirty = 1;
    end
    
    properties
        
        %STARTB List of field point coordinates, Zeeman splitting in MHz.
        StartB = [];
        
        %MESHITERATIONS Number of iterations to mesh calculation.
        MeshIterations = 3;
        
        %SORT Sort the eigenenergies by eigenvector similarity.
        Sort = 0;
        
        %QUADGRAD Quad gradient, B_Z = -2 * B' z. B' is in Gauss/cm.
        QuadGrad = 0;
        
        Verbose = 0;
        
    end
    
    methods
        
        function instance = ZAxisSampler(apCalc)
            %ZAXISSAMPLER Creates a new ZAxisSampler instance.
            %   The constructor must be given a reference to an APCalculator
            %   object which will be used to sample the field points.
            
            if ~isa(apCalc, 'AP.Calculator')
                error('apCalc must be an AP.Calculator object.');
            end
            
            instance.APCalculator = apCalc;
        end
        
        Sample(instance);
        
        function zs = Z(instance)
            %Z Get the z-axis coordinates of this z-axis sampler.
            %   Requires a quadrupole gradient to be specified.
            
            if (instance.Dirty)
                instance.dirtyError()
            end
            
            if (instance.QuadGrad <= 0)
                error('Quadrupole gradient must be > 0 for ZAxisSampler to map Zeeman splittings to spatial locations.');
            end
        end
        
        function Bs = B(instance)
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
        
        function set.StartB(instance, val)
            instance.StartB = val;
            instance.Dirty = 1;
        end
        
        function set.MeshIterations(instance, val)
            instance.MeshIterations = val;
            instance.Dirty = 1;
        end
        
        function set.QuadGrad(instance, val)
            instance.QuadGrad = val;
            instance.Dirty = 1;
        end
        
        function set.Sort(instance, val)
            instance.Sort = val;
            instance.Dirty = 1;
        end
        
        function dirtyError(~)
            error('Properties of the sampler have changed since calculation, or the eigenenergies have not been calculated yet.');
        end
    end
    
end

