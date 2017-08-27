classdef TestF1General < matlab.unittest.TestCase
    %TESTF1GENERAL Unit tests for the F1General Hamiltonian.
    
    methods (Test)
        
        function TestAgreesWithGeneral(testCase)
            %TESTAGREESWITHGENERAL
            %
            %   Tests that the F=1 Hamiltonian agrees with the
            %   general F Hamiltonian.
            
            omega0 = 1; RFs = 2; phase = 0;
            
            for theta = [ 0 0.3 1 ]
                for gamma = [ 0 0.3 1 ]
                    for gFuBBx = [ 0 0.5 1 ]
                        for gFuBBy = [ 0 0.5 1 ]
                            for gFuBBz = [ 0 0.5 1 ]
                                for py = [ 0 1 ]
                                    for pz = [ 0 1 ]
                                        approxH = @(t) AP.Hamiltonian.F1General( t, omega0, RFs, theta, gamma, gFuBBx, gFuBBy, gFuBBz, py, pz, phase );
                                        GeneralFH = @(t) AP.Hamiltonian.General( t, 1, omega0, RFs, theta, gamma, gFuBBx, gFuBBy, gFuBBz, py, pz, phase );
                                        
                                        for t=[ 1 2 3 ]
                                            testCase.assertEqual(approxH(t), GeneralFH(t), 'AbsTol', 1e-5);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
    end
    
end

