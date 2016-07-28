function [freqs, atomN] = rabispectra(exRabis, exN, Rabis, spectra)
% What I want is a way to map a given Rabi freq and atom number remaining
% to all possible rfs that rabi freq could be.

% spectra is a matrix with transition index along first index and the Rabi
% frequencies for these values along second index.
% 
% Different transitions
% ^
% |
% |
% |----------------> Different rabi freqs
% 
% Rabi freqs are specified by Rabi.
% 
% I would like to take in a list of rabi frequencies and lifetimes
% associated with them, then use spectra to map those lifetimes into
% transition frequencies. I'm not immortal and as a result I'm going to do
% a quick and dirty fix for this, because I'm likely only going to run this
% code once.

%exRabis = data42(:,1) * BaseRabi40;
%exN     = data42(:,2);

freqs = [];
atomN = [];

% enumerate over each rabi transition.
for ti=1:size(spectra, 1)
   
    % create interpolating function from spectra as LUT
    try
        f = interp1(Rabis, spectra(ti, :), exRabis);
        [~,i] = sort(f);
        freqs = [ freqs f(i)' NaN ];
        atomN = [ atomN exN(i)' NaN ];
    catch e 
        warning(e.message);
    end
    
end

end
