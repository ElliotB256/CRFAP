classdef (Abstract) AbstractSampler < handle
    %ABSTRACTSAMPLER Abstract base class for all samplers
    %   Base class for all samplers. A sampler provides a means of
    %   selecting points at which the eigenenergies should be calculated.
    %   
    
    properties
    end
    
    methods (Abstract)
        Sample(instance);
    end
    
end

