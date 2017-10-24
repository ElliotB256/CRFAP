classdef (Abstract) AbstractSampler < handle
    %ABSTRACTSAMPLER Abstract base class for all samplers
    %   Base class for all samplers. A sampler provides a means of
    %   selecting points at which the eigenenergies should be calculated.
    %
    
    properties
        
        %APCALCULATOR Reference to APCalculator used for field point calc.
        APCalculator;
        
    end
    
    methods (Abstract)
        Sample(instance);
    end
    
    methods
        
        function instance = AbstractSampler(calculator)
            
            if ~isa(calculator, 'AP.Calculator')
                error('apCalc must be an AP.Calculator object.');
            end
            
            instance.APCalculator = calculator;
        end
        
    end
    
end

