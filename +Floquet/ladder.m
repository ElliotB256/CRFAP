function [ lad ] = ladder( RFs, n, F )
%LADDER Arranges given energy states in a ladder.
% Syntax: ladder( RFs, n, F )
%  RFs: Radio frequencies of the system, used to determine system 
%       periodicity.
%  n : number of times to ladder
%  F : AP energy levels

fundamental = MRF.GetFundamental(RFs);
offsets = fundamental*(1:n)';
offsets = kron(offsets, ones(size(F,1), 1));
lad = repmat(F, n, 1) + repmat(offsets, 1, size(F, 2)); 

end

