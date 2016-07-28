HighRF = 4.2;
n = [ 7:13 20];
freqSep = repmat(HighRF, 1, length(n)) ./ n;

RFs = {};
for df=freqSep
   RFs{end+1} = [ -2*df; -df; 0] + HighRF; 
end

amps = MRF.Calculations.getFlatAmplitudes(RFs, 100);