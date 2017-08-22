function [ H ] = F1LinPolBottomOfShell( t, omega0, RFs, gFuBB, phase )
%F1LINPOLBOTTOMOFSHELL Get H(t) at bottom of shell for lin polarised RF.
% 
% gFuBB: the amplitude of the rf field expressed in MHz

% THIS FUNCTION WILL BE CALLED FREQUENTLY. Make it as lean as possible.
c = gFuBB;

H0 = [ 
      -omega0,      0,      0,       ;
            0,      0,      0,       ;
            0,      0, omega0,       ;
     ];

% Calculate coherence terms
st    = sum( c .* sin(t .* RFs + phase), 1) / (2.^0.5);

Hc = [     ...
         0,     st,      0,       ;
        st,      0,     st,       ;
         0,     st,      0,       ;
     ];
    
H = Hc + H0;

end