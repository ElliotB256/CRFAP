classdef LineSampler < AP.Sampler.AbstractLineSampler
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
    
    properties (Dependent)
        
        %B Get the field coordinates of this z-axis sampler (MHz).
        B;
       
        %POTENTIALENERGIES Eigenenergy including gravitational sag
        PotentialEnergies;
        
    end
    
    properties
        
        %STARTB List of field point coordinates, Zeeman splitting in MHz.
        StartB = [];
        
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
            
            instance = instance@AP.Sampler.AbstractLineSampler(calculator);
            instance.Theta = theta;
            instance.Gamma = gamma;
        end
        
        %Sample(instance);
        
        function coords = GetCoords(instance)
            %GETCOORDS Returns the coordinates at each field point.
            coords = instance.TransformThetaGamma2XYZ(instance.B, instance.Theta, instance.Gamma);
        end
        
        function Bs = get.B(instance)
            %B Get the field coordinates of this z-axis sampler.
            %   Returns zeeman splitting in MHz.
            
            if (instance.Dirty)
                instance.dirtyError()
            end
            
            % Convert lambdas to B 
            [Bs,~,~] = instance.lambdaToBTG(instance.mLambda);
        end
        
        function set.StartB(instance, val)
            instance.StartB = val;
            instance.Dirty = 1;
            instance.InitialLambda = (val - min(val))./(max(val)-min(val));
        end
        
        function Esag = get.PotentialEnergies(instance)
           %POTENTIALENERGIES Get eigenenergies plus gravitational sag
            Esag = instance.Eigenenergies + Util.gpe(instance.Z, instance.APCalculator.Atom.Mass);
        end

        function set.Theta(instance, val)
            instance.Theta = val;
            instance.Dirty = 1;
        end
        
        function set.Gamma(instance, val)
            instance.Gamma = val;
            instance.Dirty = 1;
        end
        
        function [B,theta,gamma] = lambdaToBTG(instance, lambda)
            B = min(instance.StartB) + (max(instance.StartB)-min(instance.StartB))*lambda;
            theta = instance.Theta;
            gamma = instance.Gamma;
        end
        
    end
    
end

