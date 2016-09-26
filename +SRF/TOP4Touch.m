function [ top ] = TOP4Touch( RF, BGrad, zsf )
%TOP4TOUCH Calculates the top field at which the resonant ellipsoids will
%touch
% RF: radiofrequency (Mhz)
% BGrad: Quadrupole strength (Gauss/cm)
% zsf: Zeeman splitting in MHz/Gauss

top = SRF.resonantEllipsoidWidth(RF, BGrad, zsf) * 1e-4 * BGrad;

end

