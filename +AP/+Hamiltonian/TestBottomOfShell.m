classdef TestBottomOfShell < matlab.unittest.TestCase
    %TESTBOTTOMOFSHELL Unit tests for general polarisation along z<0.
    
    methods (Test)
        
        function TestAgreesWithGeneral(testCase)
            %TESTAGREESWITHGENERAL
            %
            %   Tests that the approximated Hamiltonian agrees with the
            %   general Hamiltonian in the limits required for this
            %   approximation (along z<0).
            
            omega0 = 1; RFs = 2;
            for F = [ 1 2 ]
                for pz = [ 0 1 ]
                    for py = [ 0 1 ]
                        for gFuBBx = [ 0 1 ]
                            for gFuBBy = [ 0 1 ]
                                for gFuBBz = [ 0 1 ]
                                    approxH = @(t) AP.Hamiltonian.BottomOfShell( t, F, omega0, RFs, gFuBBx, gFuBBy, gFuBBz, py, pz, 0 );
                                    H = @(t) AP.Hamiltonian.General( t, F, omega0, RFs, 0, 0, gFuBBx, gFuBBy, gFuBBz, py, pz, 0 );
                                    
                                    for t=[ 1 2 3 ]
                                        testCase.assertEqual(H(t), approxH(t), 'AbsTol', 1e-5);
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
