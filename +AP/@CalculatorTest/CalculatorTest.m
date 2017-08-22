classdef CalculatorTest < matlab.unittest.TestCase
    %CALCULATORTEST Units tests for the AP Calculator class.
    
    methods (Test)
        
        function testCircPolSymmetry_ShellBottom_F1(testCase)
            %TESTCIRCPOLSYMMETRY
            %
            %   The circ dressing rf hamiltonian should
            %   be of the form [ -x . 0 ; . 0 . ; 0 . x ] for F=1.
            
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
    end
end