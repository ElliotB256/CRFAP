classdef Calculator < handle
    %CALCULATOR Calculator used for calculating properties of the TAAP.
    
    properties
        
        %ATOM Species used in the TAAP
        Atom;
        
        %BTOP TOP field amplitude, (G)
        BTop
        
        %RF Dressing amplitude (G)
        BRF;
        
        %RF Dressing frequency (MHz)
        RF;
        
        %QUADGRAD Quadrupole field gradient (G/cm)
        QuadGrad
        
        %U Potential energy in the TAAP
        U
        
        %ANISOTROPY Anisotropy of the TOP field, By=Anisotropy*Bx
        Anisotropy = 1;
        
        %TIMEAVERAGESTEPS Number of evaluations for time averaging
        TimeAverageSteps = 10;
        
    end
    
    properties (SetAccess=protected)
       
        %MFTILDE Dressed state of atoms
        MFTilde = [];
        
    end
    
    methods
        
        function context = Calculator(RF, BRF, Btop, QuadGrad)
           %CALCULATOR Create a new instance of TAAP.Calculator
           
            context.RF = RF;
            context.BRF = BRF;
            context.BTop = Btop;
            context.QuadGrad = QuadGrad;
            context.WithSpecies(87);
            
        end
        
        function context = WithSpecies(context, spec, varargin)
            
            ip = inputParser;
            ip.addParameter('mFTilde', []);
            ip.parse(varargin{:});
            
            switch (spec)
                case 85
                    F = 2;
                    gFuB = 0.7*2/3;
                    M = 85;
                case 87
                    F = 1;
                    gFuB = 0.7;
                    M = 87;
                otherwise
                    error('unknown species');
            end
            
            context.MFTilde = F;
            if ~isempty(ip.Results.mFTilde)
               context.MFTilde = ip.Results.mFTilde; 
            end
            context.Atom = AP.Atom(F, gFuB, M);
            
        end
        
        function context = WithAnisotropy(context, anis)
            
            context.Anisotropy = anis;
            
        end
        
        function pot = Calculate(context, x, y, z)
            
            % Get displacement of quadrupole centre in microns
            xTOP = context.BTop/context.QuadGrad * 1e4;
            
            % Define potential
            trap = @(a,b,c,nt) SRF.ShellTrap(...
                a - xTOP .* cos(2*pi*nt), ...
                b - context.Anisotropy * xTOP .* sin(2*pi*nt), ...
                c,...
                abs(context.Atom.gFuB), context.QuadGrad, context.RF, context.BRF, context.MFTilde) + ...
                Util.gpe(c, context.Atom.Mass);
            
            pot = Util.timeAverage(@(t) trap(x,y,z,t), context.TimeAverageSteps);
            
        end
        
    end
    
end

