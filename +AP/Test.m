%% Test
% Runs all unit tests

import matlab.unittest.TestSuite

suite = TestSuite.fromPackage('AP.Hamiltonian');
result = run(suite);