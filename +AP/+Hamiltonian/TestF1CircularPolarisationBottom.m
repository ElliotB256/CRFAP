classdef TestF1CircularPolarisationBottom < matlab.unittest.TestCase
    %TESTF1CIRCULARPOLARISATIONBOTTOM Unit tests for circular polarised dressing RF
    
    methods (Test)
        
        function TestSymmetry_ShellBottom_F1(testCase)
            %TestSymmetry_ShellBottom_F1
            %
            %   Test symmetries of the dressing RF Hamiltonian.
            
            RF = 3; BRF = 1; omega0 = 0;
            ap = AP.Calculator().CircularPolarised(RF, BRF);
            H = ap.GetHamiltonian(omega0, 0, 0);
            
            testTimes = [ 1 2 3 ];
            for t=testTimes
                h = H(1);
                % Coherences on one diagonal should equal each other, should
                % equal conjugate of coherences on the other diagonal
                testCase.verifyTrue(h(1,2) == h(2,3) && h(2,1) == h(3,2));
                testCase.verifyTrue(h(1,2) == conj(h(2,1)) && h(2,1) == conj(h(1,2)));
                
                % Test that diagonal elements are zero for omega0 = 0
                testCase.verifyTrue(h(1,1) == 0 && h(2,2) == 0 && h(3,3) == 0);
                
                % Test that second-diagonal terms are zero.
                testCase.verifyTrue(h(3,1) == 0 && h(1,3) == 0);
            end
        end
        
        function TestOffDiagonalTermsAreZeroWithNoRF(testCase)
            %TESTOFFDIAGONALTERMSAREZEROWTIHNORF
            %
            %   Off-diagonal terms should be zero with no applied rf.
            
            RF = 3; B = 0;
            ap = AP.Calculator().CircularPolarised(RF, B);
            H = ap.GetHamiltonian(1, 0, 0);
            for t = [ 0 1 2 ]
                h = H(t);
                testCase.verifyEqual(h(1,2), 0, 'RelTol', 1e-5);
                testCase.verifyEqual(h(1,3), 0, 'RelTol', 1e-5);
                testCase.verifyEqual(h(2,1), 0, 'RelTol', 1e-5);
                testCase.verifyEqual(h(2,3), 0, 'RelTol', 1e-5);
                testCase.verifyEqual(h(3,1), 0, 'RelTol', 1e-5);
                testCase.verifyEqual(h(3,2), 0, 'RelTol', 1e-5);
            end
        end
        
        function TestDiagonalTerms(testCase)
            %TESTDIAGONALTERMS
            %
            %   Diagonal terms should be equal to n*omega0*sign(gF), with
            %   n=-1...1. This comes from the diagonal terms equal to
            %   gFuBB, with B equal to omega0/gFuB.
            
            RF = 3; BRF = 1;
            ap = AP.Calculator().CircularPolarised(RF, BRF);
            
            for omega0 = [ 0 1 ]
                H = ap.GetHamiltonian(omega0, 0, 0);
                h = H(1);
                testCase.verifyEqual(omega0*sign(ap.Atom.gFuB) * (-1), h(1,1), 'AbsTol', 1e-5);
                testCase.verifyEqual(omega0*sign(ap.Atom.gFuB) * ( 0), h(2,2), 'AbsTol', 1e-5);
                testCase.verifyEqual(omega0*sign(ap.Atom.gFuB) * (+1), h(3,3), 'AbsTol', 1e-5);
            end
        end
        
        function TestAgreesWithGeneral(testCase)
            %TESTAGREESWITHGENERAL
            %
            %   Tests that the circular bottom Hamiltonian agrees with the
            %   general Hamiltonian in the limits required for this
            %   approximation.
            
            omega0 = 1; RFs = 2;
            
            for gFuBB = [ 0.5 1 ]
                approxH = @(t) AP.Hamiltonian.F1CircPolBottomOfShell( t, omega0, RFs, gFuBB, 0 );
                H = @(t) AP.Hamiltonian.F1General( t, omega0, RFs, 0, 0, gFuBB, gFuBB, 0, pi/2, 0, 0 );
                
                for t=[ 1 2 3 ]
                    testCase.verifyEqual(H(t), approxH(t), 'AbsTol', 1e-5);
                end
            end
        end
        
    end
    
end

