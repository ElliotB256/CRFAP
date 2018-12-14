classdef (Abstract) QuadrupoleSampler < AP.Sampler.AbstractSampler & matlab.mixin.SetGet
    %QUADRUPOLESAMPLER Samples dressed potentials in a static quad field.
    %   The quadrupole sampler provides mapping functions to transform from
    %   x,y,z to rotation angles theta, gamma for the local quantisation
    %   axis.
    
    properties (SetAccess=protected, GetAccess=protected)
       
        %FILEDIRTY Represents the transient property 'Dirty' for save/load.
        FileDirty;
        
    end
    
    properties (Transient,SetAccess=protected)
        
        %DIRTY The instance is dirty if properties are changed/before calc.
        Dirty = 1;
        
    end
    
    properties (Dependent)
       
        %X Get the y-axis coordinates of this line sampler (microns).
        X;
        
        %Y Get the y-axis coordinates of this line sampler (microns).
        Y;
        
        %Z Get the y-axis coordinates of this line sampler (microns).
        Z;
        
    end
    
    properties
       
        %QUADGRAD Quad gradient, B_Z = -2 * B' z. B' is in Gauss/cm.
        QuadGrad = 0;
        
    end
    
    methods
        
        function instance = QuadrupoleSampler(calculator)
           instance = instance@AP.Sampler.AbstractSampler(calculator);
        end
        
        function set.QuadGrad(instance, val)
            instance.QuadGrad = val;
            instance.Dirty = 1;
        end
        
        function dirtyError(~)
            error('Properties of the sampler have changed since calculation, or the eigenenergies have not been calculated yet.');
        end
        
        function coords = TransformThetaGamma2XYZ(instance, zeemanSplit, theta, gamma)
            %GETCOORDS Returns the coordinates at each field point.
            %   zeemanSplit is the zeeman splitting in MHz. theta and gamma
            %   are the rotation angles for the quantisation axis.
            
            if (instance.QuadGrad <= 0)
                error('Quadrupole gradient must be > 0 for ZAxisSampler to map Zeeman splittings to spatial locations.');
            end
            
            % get uBgF in MHz/Gauss
            gFuB = abs(instance.APCalculator.Atom.gFuB);
            
            
            % The quantisation axis fixes the ratio x:y:z. The magnitude of
            % the zeeman splitting determines the position along this ray.
            mag = zeemanSplit/(((instance.QuadGrad))*gFuB)*1e4; % microns
            
            z = mag .* cos(theta) / 2;
            x = mag .* sin(theta) .* cos(gamma);
            y = mag .* sin(theta) .* sin(gamma);
            
            coords = struct('x', x(:)', 'y', y(:)', 'z', z(:)');
            
        end
        
        function xs = get.X(instance)
            %X Get the x-axis coordinates of this line sampler (microns).
            %   Requires a quadrupole gradient to be specified.
            coords = instance.GetCoords();
            xs = coords.x;
        end
        
        function ys = get.Y(instance)
            %Y Get the y-axis coordinates of this line sampler (microns).
            %   Requires a quadrupole gradient to be specified.
            coords = instance.GetCoords();
            ys = coords.y;
        end
        
        function zs = get.Z(instance)
            %Z Get the z-axis coordinates of this z-axis sampler (microns).
            %   Requires a quadrupole gradient to be specified.
            coords = instance.GetCoords();
            zs = coords.z;
        end
        
    end
    
    methods (Static)
       
        function b = loadobj(a)
            for i=1:length(a)
                b(i) = a(i);
                b(i).Dirty = a(i).FileDirty; 
            end
        end
        
        function b = saveobj(a)
            for i=1:length(a)
                b(i) = a(i);
                b(i).FileDirty = a(i).Dirty;
            end
        end
        
    end
    
end

