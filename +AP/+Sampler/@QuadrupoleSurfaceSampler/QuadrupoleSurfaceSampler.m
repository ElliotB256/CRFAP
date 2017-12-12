classdef QuadrupoleSurfaceSampler < AP.Sampler.QuadrupoleSampler
    %QUADRUPOLESURFACESAMPLER Samplers points on resonant spheroid
    
    properties
        
        %ALPHA Vertical angles alpha at which to sample
        Alpha;
        
        %BETA Radial angles at which to sample
        Beta;
        
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
        
        %GAMMA Rotation angles gamma for frame transformation;
        Gamma;
        
        %THETA Rotation angles theta for frame transformation;
        Theta;
        
        %EIGENVECTORS Rank 1,2 are eigenvectors, size(EV,3) == length(B)
        Eigenvectors
        
        %ENERGIES Get the energies of this z-axis sampler at field points.
        Eigenenergies;
        
        %POTENTIALENERGIES Eigenenergy including gravitational sag
        PotentialEnergies;
        
    end
    
    methods
        
        function instance = QuadrupoleSurfaceSampler(calculator, alpha, beta)
            %QUADRUPOLESURFACESAMPLER Creates a new QuadrupoleSurfaceSampler instance.
            %   The constructor must be given a reference to an APCalculator
            %   object which will be used to sample the field points.
            
            if nargin < 3
                beta = linspace(0, 2 * pi, 12);
                beta = beta(1:end-1);
            end
            
            if nargin < 2
                alpha = linspace(0, pi/2, 6);
            end
            
            instance = instance@AP.Sampler.QuadrupoleSampler(calculator);
            instance.Alpha = alpha;
            instance.Beta = beta;
            instance.RF = calculator.RF(1);
        end
        
        function alpha = GetAlphas(instance)
            a = instance.Alpha(:);
            b = instance.Beta(:);
            alpha = repmat(a', length(b), 1); alpha = alpha(:);
            clear a b
        end
        
        function beta = GetBetas(instance)
            a = instance.Alpha(:);
            b = instance.Beta(:);
            beta = repmat(b, 1, length(a)); beta = beta(:);
            clear a b
        end
        
%         function coords = GetCoords(instance)
%             %GETCOORDS Returns the coordinates at each field point.
%             theta = instance.Theta;
%             gamma = instance.Gamma;
%             coords = cell(length(gamma),1);
%             for i=1:length(theta)
%                 coords{i} = instance.TransformThetaGamma2XYZ(instance.RF, theta(i), gamma(i));
%             end
%             c = cat(1,coords{:});
%             coords = struct(...
%                 'x', cat(1, c.x), ...
%                 'y', cat(1, c.y), ...
%                 'z', cat(1, c.z) ...
%                 );
%         end
        

function c = GetNCoords(instance)
    c = GetCoords(instance);
end
% prev NCoords
        function coords = GetCoords(instance)
            
            alpha = instance.GetAlphas();
            beta = instance.GetBetas();
            
            x = zeros(1, length(alpha));
            y = zeros(1, length(alpha));
            z = zeros(1, length(alpha));
%             for i=1:length(alpha)
                x = sin(alpha) .* cos(beta);
                y = sin(alpha) .* sin(beta);
                z = cos(alpha) / 2;
%             end
            
            coords = struct('x', x, 'y', y, 'z', z);
            
        end
        
        function gamma = get.Gamma(instance)
            coords = instance.GetNCoords();
            gamma = acos( (coords.x.^2 + 4*coords.z.^2).^0.5 ./ (coords.x.^2 + coords.y.^2 + 4 * coords.z.^2).^0.5 );
        end
        
        function theta = get.Theta(instance)
            coords = instance.GetNCoords();
            theta = acos( (-2.*coords.z) ./ (coords.x.^2 + 4 * coords.z.^2).^0.5 );
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
        
        Sample(instance);
        
    end
    
end

