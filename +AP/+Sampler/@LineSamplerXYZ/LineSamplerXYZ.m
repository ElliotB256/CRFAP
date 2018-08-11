classdef LineSamplerXYZ < AP.Sampler.AbstractLineSampler
    %LINESAMPLERXYZ Samples the adiabatic potential along a line.
    %   Calculates eigenenergies and vectors of the dressed potential along
    %   a line between InitialSamplePoints(1) and (end). Employs meshing to
    %   refine the calculation, up to a specified number of iterations.
    %
    %   A quadrupole gradient can be specified, in which case the meshing
    %   will also account for gravitational sag.
    
    properties (Dependent)
        
    end
    
    properties
        
        %INITIALSAMPLEPOINTS Initial points for sampling
        InitialSamplePoints;
        
    end
    
    methods
        
        function instance = LineSamplerXYZ(calculator)
            %LINESAMPLER3D Creates a new LineSampler3D instance.
            %   The constructor must be given a reference to an APCalculator
            %   object which will be used to sample the field points.
            instance = instance@AP.Sampler.AbstractLineSampler(calculator);
        end
        
        function set.InitialSamplePoints(instance, val)
            
            if size(val, 2) ~= 3
                error('InitialSamplePoints must be specified as rows of X,Y,Z.');
            end
            if size(val, 1) < 2
                error('At least two points must be specified.');
            end
            
            instance.InitialSamplePoints = val;
            instance.Dirty = 1;
            
            x = instance.InitialSamplePoints(:,1);
            y = instance.InitialSamplePoints(:,2);
            z = instance.InitialSamplePoints(:,3);
            
            % Initial sample points may not be evenly spaced. The variable, lambda, is
            % equivalent to displacement along the line.
            cx = max(x) - min(x);
            cy = max(y) - min(y);
            cz = max(z) - min(z);
            if cx > cy && cx > cz
                l = x;
            elseif cy > cz
                l = y;
            else
                l = z;
            end
            lambda = (l - l(1))./(l(end)-l(1));
            instance.InitialLambda = lambda';
        end
        
        function instance = withLinearSpan(instance, startP, endP, N)
            %WITHLINEARSPAN Sample N points initially from startP to endP.
            %  The N points are linearly spaced along the line between the
            %  points startP and endP.
            
            if ~all(size(startP) == [ 1 3 ])
                error('startP must be a row of [ x y z ].');
            end
            
            if ~all(size(endP) == [ 1 3 ])
                error('endP must be a row of [ x y z ].');
            end
            
            x = linspace(startP(1), endP(1), N)';
            y = linspace(startP(2), endP(2), N)';
            z = linspace(startP(3), endP(3), N)';
            
            instance.InitialSamplePoints = [ x y z ];
            
        end
        
        function [B,theta,gamma] = lambdaToBTG(instance, lambda)
            
            x0 = instance.InitialSamplePoints(1,1);
            y0 = instance.InitialSamplePoints(1,2);
            z0 = instance.InitialSamplePoints(1,3);
            
            dx = instance.InitialSamplePoints(end,1) - instance.InitialSamplePoints(1,1);
            dy = instance.InitialSamplePoints(end,2) - instance.InitialSamplePoints(1,2);
            dz = instance.InitialSamplePoints(end,3) - instance.InitialSamplePoints(1,3);
            
            x = x0 + dx*lambda;
            y = y0 + dy*lambda;
            z = z0 + dz*lambda;
            
            B = abs(instance.APCalculator.Atom.gFuB)*(instance.QuadGrad * 1e-4) .* (x.^2 + y.^2 + 4 * z.^2).^0.5;
            %gamma = acos(x ./ (x.^2 + y.^2).^0.5);
            gamma = atan2(y, x);
            theta = acos(2.*z ./ (x.^2 + y.^2 + 4 * z.^2).^0.5);
            
            gamma(isnan(gamma)) = 0; % z-axis
            theta(isnan(theta)) = 0; % origin
            
        end
        
        function h = IsHorizontal(instance)
            %ISHORIZONTAL Returns true if line is horizontal
            h = instance.InitialSamplePoints(1,3) == instance.InitialSamplePoints(end,3);
        end
        
    end
    
end