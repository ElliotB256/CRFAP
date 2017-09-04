classdef TestF1CircularPolarisation < matlab.unittest.TestCase
    %TESTF1CIRCULARPOLARISATION Unit tests for circular polarised dressing RF
    
    methods (Test)
        
        function TestAgreesWithGeneral(testCase)
            %TESTAGREESWITHGENERAL
            %
            %   Tests that the circular bottom Hamiltonian agrees with the
            %   general Hamiltonian in the limits required for this
            %   approximation.
            
            omega0 = 1; RFs = 2;
            
            for theta = [ 0 0.5 1.0 1.5 pi ]
                for gFuBB = [ 0.5 1 ]
                    approxH = @(t) AP.Hamiltonian.F1CircPol( t, omega0, RFs, gFuBB, theta, 0 );
                    H = @(t) AP.Hamiltonian.F1General( t, omega0, RFs, theta, 0, gFuBB, gFuBB, 0, pi/2, 0, 0 );
                    
                    for t=[ 1 2 3 ]
                        testCase.assertEqual(H(t), approxH(t), 'AbsTol', 1e-5);
                    end
                end
            end
        end
        
    end
    
end

