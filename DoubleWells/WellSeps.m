

n = [ 5 6 7 ]';
df = 0.6;
RFs = n*df;
Rabi = df * [ 0.7 0.3 0.7 ]';
qdrpGrad = 200;

CalcWellSeparations

% Small df scheme
n = [ 40 41 42 ]';
df = 0.1;
RFs = n*df;
Rabi = df * [ 0.7 0.3 0.7 ]';
qdrpGrad = 200;

CalcWellSeparations

% Really push qdrp
n = [ 40 41 42 ]';
df = 0.1;
RFs = n*df;
Rabi = df * [ 0.7 0.3 0.7 ]';
qdrpGrad = 400;

CalcWellSeparations