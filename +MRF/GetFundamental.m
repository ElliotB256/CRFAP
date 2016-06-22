function [ freq ] = GetFundamental( RFs )
%GETFundamental Calculates periodicity of given RF frequencies. Assumes
%frequencies have common difference. Only works for RFs with a common
%fundamental.

if length(RFs) == 1
    freq = RFs;
else
    
    RFs = sort(RFs);
    if any(diff(diff(RFs)))
        error('RFs must have common difference');
    end
    
    d = RFs(2) - RFs(1);
    
    % divide RFs by d.
    n = round(RFs / d);
    
    test = n*d-RFs;
    if any(abs(test) > 1e-10)
        error('RFs must have common fundamental');
    end
    
    % get gcd of numbers
    % b = MRF.lcms(n)
    % for now: assume gcd is always just '1'. This is only a performance
    % concern.
    freq = d;
    
    
    % t = (0:0.01/freq:1/freq) * 2 * pi;
    % plot(t, [cos(t * RFs(1)); cos(t * RFs(2)); cos(t * RFs(3))]);
end

end

