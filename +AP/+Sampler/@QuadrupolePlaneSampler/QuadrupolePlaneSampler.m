classdef QuadrupolePlaneSampler < AP.Sampler.QuadrupoleSampler
    %QUADRUPOLEPLANESAMPLER Calculate field points in the xz-plane
    %   Calculates field points for the quadrupole field that are in the
    %   plane y=0. This sampler employs a fan of line samplers to calculate
    %   the dressed states, and as such can make use of meshing to improve
    %   accuracy.
    %
    %   This sampler uses a static quadrupole field for the static field
    %   components.
    
    properties (SetAccess=protected)
       
        %LINESAMPLERS LineSampler objects used by the plane sampler to
        %evaluate potentials.
        LineSamplers = [];
        
    end
    
    properties
       
        %STARTB List of starting coordinates for rays, Zeeman splitting in MHz.
        StartB = [];
        
        %MESHITERATIONS Number of iterations to mesh calculation.
        MeshIterations = 3;
        
        %SORTLINES Sort the eigenenergies in rays by eigenvector similarity.
        SortLines = 1;
        
        %VERBOSE Should the sampler print debug information?
        Verbose = 1;
        
        %GAMMA Constant gamma used for the plane.
        Gamma;
        
        %THETARANGE Range of theta values to evaluate potential over.
        ThetaRange; % must be length(2), 0 to pi
        
        %RAYNUMBER Number of rays to cast for potential calculation.
        RayNumber; % must be integer > 0
    end
    
    methods
        
        function instance = QuadrupolePlaneSampler(calculator, thetaRange, rayNumber, gamma)
            %QUADRUPOLEPLANESAMPLER Creates a new QuadrupolePlaneSampler instance.
            %   The constructor must be given a reference to an APCalculator
            %   object which will be used to sample the field points.
            
            if nargin < 4
                gamma = 0;
            end
            
            if nargin < 3
                rayNumber = 10;
            end
            
            if nargin < 2
                thetaRange = [ 0 pi ];
            end
            
            instance = instance@AP.Sampler.QuadrupoleSampler(calculator);
            instance.RayNumber = rayNumber;
            instance.ThetaRange = thetaRange;
            instance.Gamma = gamma;
        end
        
        function coords = GetCoords(instance)
            %GETCOORDS Returns the coordinates at each field point.
            
            if (instance.Dirty)
                instance.dirtyError();
            end
            
            % Get coordinates of each line sampler
            s = struct('x', {}, 'y', {}, 'z', {});
            for i=1:length(instance.LineSamplers)
               ls = instance.LineSamplers{i};
               s(end+1) = ls.GetCoords();
            end
            
            % Bash coords into list of coords
            coords = struct('x', [s.x], 'y', [s.y], 'z', [s.z]);
        end
        
        function xs = X(instance)
            %X Get the x-axis coordinates of this plane sampler (microns).
            %   Requires a quadrupole gradient to be specified.
            coords = instance.GetCoords();
            xs = coords.x;
        end
        
        function ys = Y(instance)
            %Y Get the y-axis coordinates of this plane sampler (microns).
            %   Requires a quadrupole gradient to be specified.
            coords = instance.GetCoords();
            ys = coords.y;
        end
        
        function zs = Z(instance)
            %Z Get the z-axis coordinates of this plane sampler (microns).
            %   Requires a quadrupole gradient to be specified.
            coords = instance.GetCoords();
            zs = coords.z;
        end
        
        function Es = E(instance)
            %E Get the energies of this plane sampler at field points.
            
            if (instance.Dirty)
                instance.dirtyError()
            end
            
            % Get energies of each line sampler
            Es = {};
            for i=1:length(instance.LineSamplers)
               ls = instance.LineSamplers{i};
               Es{end+1} = ls.E();
            end
            
            % Concatenate into single array of energy
            Es = [Es{:}];
        end
        
    end
    
end

