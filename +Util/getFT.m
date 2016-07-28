function [ frequencies, amplitudes ] = getFT( time, quantity )
%GETFT Performs a fourier transform of the given quantity.
% Example:
% [frequencies, amplitudes] = getFT(times, quantity);

if length(time) ~= length(quantity)
   error('time and quantity must be same length!'); 
end

timestep = time(2)-time(1);
fs = 1/timestep;

%produce FFT of the data.

%Step 1 - kill the DC.
quantity = quantity - mean(quantity);

%number of bins
m = length(quantity);
n = pow2(nextpow2(m));
frequencies = (0:n-1)*(fs/n);

%fast ft for the data
amplitudes = fft(quantity,n)/length(quantity)*2;

%clip
i=1:floor(n/2);
amplitudes = amplitudes(i);
frequencies = frequencies(i);

end

