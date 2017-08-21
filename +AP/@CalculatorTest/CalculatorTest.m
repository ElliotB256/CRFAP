classdef CalculatorTest < matlab.unittest.TestCase
    %CALCULATORTEST Units tests for the AP Calculator class.
    
    methods (Test)
        
        function testLinearSymmetry_ShellBottom_F1(testCase)      
            %TESTLINEARSYMMETRY 
            % 
            %   The linear dressing rf hamiltonian should
            %   be of the form [ -x . 0 ; . 0 . ; 0 . x ] for F=1.
            
            RF = 3; BRF = 1; omega0 = 0;
            ap = AP.Calculator().LinearPolarised(RF, BRF);
            H = ap.GetHamiltonian(omega0, 0, 0);
            h = H(1);
            
            % Test that coherences are equal
            testCase.verifyTrue(h(1,2) == h(2,1) && h(1,2) == h(3,2) && h(1,2) == h(2,3));
            
            % Test that diagonal elements are zero for omega0 = 0
            testCase.verifyTrue(h(1,1) == 0 && h(2,2) == 0 && h(3,3) == 0);
            
            % Test that most diagonal terms are zero.
            testCase.verifyTrue(h(3,1) == 0 && h(1,3) == 0);
        end
    end
end