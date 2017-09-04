classdef LineSampler < AP.Sampler.AbstractSampler & matlab.mixin.SetGet
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
        Sort = 1;
        
        %QUADGRAD Quad gradient, B_Z = -2 * B' z. B' is in Gauss/cm.
        QuadGrad = 0;
        
        %VERBOSE Should the line sampler print debug information?
        Verbose = 1;
        
        %THETA Constant theta that describes this line.
        Theta = 0;
        
        %GAMMA Constant gamma that describes this line.
        Gamma = 0;
        
    end
    
    methods
        
        function instance = LineSampler(apCalc, theta, gamma)
            %LINESAMPLER Creates a new LineSampler instance.
            %   The constructor must be given a reference to an APCalculator
            %   object which will be used to sample the field points.
            
            if nargin < 3
                gamma = 0;
            end
            
            if nargin < 2
                theta = 0;
            end
            
            if ~isa(apCalc, 'AP.Calculator')
                error('apCalc must be an AP.Calculator object.');
            end
            
            instance.APCalculator = apCalc;
            instance.Theta = theta;
            instance.Gamma = gamma;
        end
        
        Sample(instance);
        
        function coords = GetCoords(instance)
            %GETCOORDS Returns the coordinates at each field point.
            %   zeemanSplitting is the zeeman splitting in MHz
            
            if (instance.Dirty)
                instance.dirtyError();
            end
            
            if (instance.QuadGrad <= 0)
                error('Quadrupole gradient must be > 0 for ZAxisSampler to map Zeeman splittings to spatial locations.');
            end
            
            % get uBgF in MHz/Gauss
            gFuB = abs(instance.APCalculator.Atom.gFuB);
            
            % rotate unit vector along z axis into a given direction.
            % Note: checked matrix has same sign as mathematica with
            % RotationMatrix[-\theta, ey];
            r1 = roty(rad2deg(instance.Theta));
            % Note: checked matrix has same sign as mathematica with
            % RotationMatrix[-\gamma, ex];
            r2 = rotx(rad2deg(instance.Gamma));
            
            rotMat = r2 * r1;
            
            % Rotate the z-axis vector to get local quantisation axis.
            lQA = rotMat * [ 0; 0; 1 ];
            
            % The quantisation axis fixes the ratio x:y:z. The magnitude of
            % the zeeman splitting determines the position along this ray.
            mag = instance.mB/(((instance.QuadGrad))*gFuB)*1e4; % microns
            a = mag * cos(instance.Gamma);
            z = -a .* cos(instance.Theta) / 2;
            x = sign(sin(instance.Theta).*cos(instance.Gamma)) * (a.^2 - 4 .* z .^2).^0.5;
            y = sign(sin(instance.Gamma)).*(mag.^2 - x.^2 - 4 * z.^2).^0.5;
            
            coords = struct('x', x, 'y', y, 'z', z);
            
        end
        
        function xs = X(instance)
            %X Get the x-axis coordinates of this line sampler (microns).
            %   Requires a quadrupole gradient to be specified.
            coords = instance.GetCoords();
            xs = coords.x;
        end
        
        function ys = Y(instance)
            %Y Get the y-axis coordinates of this line sampler (microns).
            %   Requires a quadrupole gradient to be specified.
            coords = instance.GetCoords();
            ys = coords.y;
        end
        
        function zs = Z(instance)
            %Z Get the z-axis coordinates of this z-axis sampler (microns).
            %   Requires a quadrupole gradient to be specified.
            coords = instance.GetCoords();
            zs = coords.z;
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
        
        function set.Theta(instance, val)
            instance.Theta = val;
            instance.Dirty = 1;
        end
        
        function set.Gamma(instance, val)
            instance.Gamma = val;
            instance.Dirty = 1;
        end
        
        function dirtyError(~)
            error('Properties of the sampler have changed since calculation, or the eigenenergies have not been calculated yet.');
        end
    end
    
end

