classdef (Abstract) QuadrupoleSampler < AP.Sampler.AbstractSampler & matlab.mixin.SetGet
    %QUADRUPOLESAMPLER Samples dressed potentials in a static quad field.
    %   The quadrupole sampler provides mapping functions to transform from
    %   x,y,z to rotation angles theta, gamma for the local quantisation
    %   axis.
    
    properties (Transient)
        
        %DIRTY The instance is dirty if properties are changed/before calc.
        Dirty = 1;
        
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
            
            % rotate unit vector along z axis into a given direction.
            % Note: checked matrix has same sign as mathematica with
            % RotationMatrix[-\theta, ey];
            r1 = roty(rad2deg(theta));
            % Note: checked matrix has same sign as mathematica with
            % RotationMatrix[-\gamma, ex];
            r2 = rotx(rad2deg(gamma));
            
            rotMat = r2 * r1;
            
            % Rotate the z-axis vector to get local quantisation axis.
            lQA = rotMat * [ 0; 0; 1 ];
            
            % The quantisation axis fixes the ratio x:y:z. The magnitude of
            % the zeeman splitting determines the position along this ray.
            mag = zeemanSplit/(((instance.QuadGrad))*gFuB)*1e4; % microns
            a = mag * cos(gamma);
            z = -a .* cos(theta) / 2;
            x = sign(sin(theta).*cos(gamma)) * (a.^2 - 4 .* z .^2).^0.5;
            y = sign(sin(gamma)).*(mag.^2 - x.^2 - 4 * z.^2).^0.5;
            
            coords = struct('x', x, 'y', y, 'z', z);
            
        end
    end
    
end

