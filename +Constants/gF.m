function [ zsf ] = gF( F, species )
% Zeeman splitting in MHz/Gauss, includes bohr magneton

if nargin < 2
    species = 'Rb87';
end

if nargin < 1
    F = 1;
end

switch species
    case 'Rb87'
      switch F
          case 1
              zsf = 0.7;
          case 2
              zsf = -0.7;
          otherwise
              error('unknown');
      end
    otherwise
        error('unknown');
end
        
        

end

