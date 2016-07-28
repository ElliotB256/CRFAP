function amps = getFlatAmplitudes(RFs, quad)
%GETFLATAMPLITUDES Calculates the Rabi frequencies for which wells are
%balanced and the barrier is just low enough that we have one minima, not
%two. It's probably not exact but it should be close enough.
%
% RFs should be a cell array of three RF combinations to find amplitudes
% for.
%
% quad should be the quadrupole gradient in G/cm

amps = {};

for i=1:length(RFs)
    fprintf('Set %d of %d:\n', i, length(RFs))
    RF = RFs{i};
    
    Bmin = min(RF(:))*0.9;
    Bmax = max(RF(:))*1.1;
    Bs = Bmin:((Bmax-Bmin)/8):Bmax;
    
    Rabi = repmat(RF(2)-RF(1), 3, 1)/4;
    
    try
        
    % This will largely be an iterative process. We start with some initial
    % Rabi frequencies equal to 1/4 the frequency spacing. This should be a
    % double well potential.
    debug = 1;
    
    balanced = 0;
    while balanced == 0
        
        [E, B] = getTrappedLevel(Bs, RF, Rabi, quad);
        
        if debug
            plot(B, E); pause(0.1);
        end
        
        % get minima
        [~,minima] = MRF.Calculations.getMinima(B, E);
        if length(minima) == 2
            % balance double well heights
            dblW = [];
            dblW(1) = minima(1);
            dblW(2) = minima(2);
            
            % adjust ratio of two Rabi freqs and try again
            wellDiff = diff(dblW);
            
            if abs(wellDiff) < 0.05
                % fine: balanced, continue;
                balanced = 1;
            else
                if wellDiff < 0
                    
                    Rabi(1) = Rabi(1) * 0.95;
                    Rabi(3) = Rabi(3) * 1.05;
                    fprintf('Raise f3...\n')
                else
                    
                    Rabi(1) = Rabi(1) * 1.05;
                    Rabi(3) = Rabi(3) * 0.95;
                    fprintf('Lower f3\n')
                end
            end
            
        else
            error('Found more/less than 2 minima!');
        end
        
    end
    
    % So we've found the values that approximately balance the two outer
    % freqs. Let's now try and increase the barrier height until we only
    % have one well in the trapped manifold.
    dblWell = 1;
    tries = 0;
    sCheck = 0;
    while dblWell
        
        [E, B] = getTrappedLevel(Bs, RF, Rabi, quad);
        
        if debug
            plot(B, E); pause(0.001);
        end
        
        [~,minima] = MRF.Calculations.getMinima(B, E);
        
        switch length(minima)
            case 1
                if sCheck
                    dblWell = 0;
                    fprintf('Found two consecutive Rabi freqs with single well.\n')
                else
                    fRabi = Rabi;
                    sCheck = 1;
                    fprintf('Found a single well, checking it is not spurious.\n')
                end
            case 2
                dblWell = 1;
                Rabi(2) = Rabi(2) * 1.1;
                
                if sCheck
                    fprintf('Prev single well was spurious.\n')
                end
                
                sCheck = 0;
            otherwise
                Rabi(2) = Rabi(2) * 1.1;
                tries = tries + 1;
                if tries > 3
                    error('Unhandled number of minima during barrier height increase.');
                end
        end
        
    end
        amps{end+1} = fRabi;
    catch e
        warning(e.message);
       amps{end+1} = []; 
    end
end