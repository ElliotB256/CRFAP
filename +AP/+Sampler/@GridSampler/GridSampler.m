classdef GridSampler < AP.Sampler.QuadrupoleSampler
    %GRIDSAMPLER Calculate eigenenergies over grid of points
    %   Calculates eigenenergies over a grid of points. The grid does not
    %   have to be uniformly spaced. It is specified by the vectors GX, GY,
    %   GZ.
    
    properties
        
        %VERBOSE Should the sampler output information as it calculates.
        Verbose = 1;
        
    end
    
    properties (SetAccess=protected)
        
        %SX Sample x points
        SX
        
        %SY Sample y points
        SY
        
        %SZ Sample z points
        SZ
        
        %B Magnitude of static field (at sampled points)
        B;
        
        %EIGENVECTORS Rank 1,2 are eigenvectors, size(EV,3) == length(B)
        Eigenvectors;
        
        %ENERGIES Get the energies of this sampler at field points.
        Eigenenergies;
        
    end
    
    properties (Dependent)
       
        %POTENTIALENERGIES Eigenenergy including gravitational sag
        PotentialEnergies;
        
    end
    
    properties (SetAccess=protected, GetAccess=protected)
        
        %THETA Theta (rotation axis of sampled points)
        Theta;
        
        %GAMMA Gamma (rotation axis of sampled points)
        Gamma;
        
    end
    
    methods
        
        Sample(sampler);
        
        function instance = GridSampler(ap, X,Y,Z)
            %GRIDSAMPLER Construct a new instance
            % X,Y,Z are the positions at which to sample
            
            instance = instance@AP.Sampler.QuadrupoleSampler(ap);
            instance.SX = X;
            instance.SY = Y;
            instance.SZ = Z;
            
        end
        
        function c = GetCoords(sampler)
            c = struct('x', sampler.SX(:), 'y', sampler.SY(:), 'z', sampler.SZ(:));
        end
        
        function E = get.PotentialEnergies(sampler)
            gpe = Util.gpe(sampler.Z', sampler.APCalculator.Atom.Mass);
            E = repmat(gpe, size(sampler.Eigenenergies, 1), 1) + sampler.Eigenenergies;
        end
    end
end

