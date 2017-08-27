classdef TestLinearPolarisation < matlab.unittest.TestCase
    %TESTLINEARPOLARISATION Unit tests for linear polarised dressing RF.
    %   These tests exploit symmetries of the Hamiltonian to test that the
    %   implementation is correct.
    
    methods (Test)
        
        function TestNoCoherenceWhenParallel(testCase)
            %TESTNOCOHERENCEWHENPARALLEL
            %   Checks that there is no coherence term when the dressing RF
            %   is parallel to the local field.
            
            % Linear field is along x. So rotate local field to be
            % parallel.
            RF = 3; B = 1;
            ap = AP.Calculator().LinearPolarised(RF, B);
            H = ap.GetHamiltonian(1, pi/2, 0);
            for t = [ 0 1 2 ]
                h = H(t);
                testCase.verifyEqual(h(1,2), 0, 'AbsTol', 1e-5);
                testCase.verifyEqual(h(1,3), 0, 'AbsTol', 1e-5);
                testCase.verifyEqual(h(2,1), 0, 'AbsTol', 1e-5);
                testCase.verifyEqual(h(2,3), 0, 'AbsTol', 1e-5);
                testCase.verifyEqual(h(3,1), 0, 'AbsTol', 1e-5);
                testCase.verifyEqual(h(3,2), 0, 'AbsTol', 1e-5);
            end
            
        end
        
        function TestCoherenceWhenPerpendicular(testCase)
            %TESTCOHERENCEWHENPERPENDICULAR
            %   Checks that coherence terms are the same magnitude when
            %   local field is rotated to be perpendicular to dressing
            %   field.
            
            % Linear field is along x. So rotate local field to be
            % parallel.
            RF = 3; B = 1;
            ap = AP.Calculator().LinearPolarised(RF, B);
            H = ap.GetHamiltonian(1, 0, 0); h1 = H(1);
            H = ap.GetHamiltonian(1, 0, pi/2); h2 = H(1);
            
            testCase.verifyEqual(h1(1,2), h2(1,2), 'AbsTol', 1e-5);
            testCase.verifyEqual(h1(1,3), h2(1,3), 'AbsTol', 1e-5);
            testCase.verifyEqual(h1(2,1), h2(2,1), 'AbsTol', 1e-5);
            testCase.verifyEqual(h1(2,3), h2(2,3), 'AbsTol', 1e-5);
            testCase.verifyEqual(h1(3,1), h2(3,1), 'AbsTol', 1e-5);
            testCase.verifyEqual(h1(3,2), h2(3,2), 'AbsTol', 1e-5);
            
        end
        
        % And all from TestLinearPolarisationBottom? Check per-Hamiltonian?
        % Check that each set of approximations agrees with the general
        % Hamiltonian in the desired limit.
        function TestAgreesWithGeneral(testCase)
            %TESTAGREESWITHGENERAL
            %
            %   Tests that the linear Hamiltonian agrees with the general
            %   Hamiltonian in the limits required for this approximation.
            
            omega0 = 1; RFs = 2;
            
            for theta = [ 0 0.3 1 ]
                for gamma = [ 0 0.3 1 ]
                    for gFuBB = [ 0.5 1 ]
                        approxH = @(t) AP.Hamiltonian.F1LinPol( t, omega0, RFs, theta, gamma, gFuBB, 0 );
                        H = @(t) AP.Hamiltonian.F1General( t, omega0, RFs, theta, gamma, gFuBB, 0, 0, 0, 0 );
                        
                        for t=[ 1 2 3 ]
                            testCase.verifyEqual(H(t), approxH(t), 'AbsTol', 1e-5);
                        end
                    end
                end
            end
        end
    end
end
