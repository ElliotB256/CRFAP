function [ gp ] = gpe( ZeemanSplit, qdrpGrad )
%GPE Calculates gravitational potential energy in MHz for given zeeman
%splitting and qdrp gradient.

BGauss = ZeemanSplit / Constants.gF;
microns = -BGauss / qdrpGrad * 1e4;
gp = gpe(microns, 87);

end

