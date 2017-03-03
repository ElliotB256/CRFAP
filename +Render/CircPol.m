function CircPol( filename, varargin )
%CIRCPOL Renders a circularly polarised rf potential in space.
% Supports multiple RFs

p = inputParser();
addRequired(p, 'filename', @ischar);
addParameter(p, 'quadGrad', 100, @isnumeric);
addParameter(p, 'rabi', 0.4, @isnumeric);
addParameter(p, 'rf', 4.2, @isnumeric);
addParameter(p, 'showGraphs', 1, @isnumeric);
addParameter(p, 'iterations', 3, @isnumeric);
addParameter(p, 'levels', 2, @isnumeric);
addParameter(p, 'gF', 0.7, @isnumeric);
addParameter(p, 'mass', 87, @isnumeric);
addParameter(p, 'F', 1, @isnumeric);
addParameter(p, 'thetas', 0:0.05:(pi/2), @isnumeric); % angles we integrate
addParameter(p, 'Bs', nan, @isnumeric); % B fields to evaluate potential at before finer meshing occurs.

p.parse(filename, varargin{:});

RF = p.Results.rf;
Rabi = p.Results.rabi;
Bs = p.Results.Bs;
if isnan(Bs)
Bs = RF-1:0.2:RF+1;
end

qdrpGrad = p.Results.quadGrad;


thetas = p.Results.thetas;

% Generate the dressing rf potential
xs = cell(1, length(thetas));
zs = cell(1, length(thetas));
vals = cell(1, length(thetas));

for i=1:length(thetas)
    theta = thetas(i);
    
    [ F, B ] = MRF.MeshedQuasiEnergies(Bs, RF, Rabi, 'iterations', p.Results.iterations, 'qdrpGrad', p.Results.quadGrad, 'F', p.Results.F, 'theta', theta, 'gF', p.Results.gF, 'mass', p.Results.mass);
    F2 = MRF.sortEnergies(B, MRF.ladder(RF, 30, F), 'F', p.Results.F);
    
    % select a level; add these values to array
    levels = p.Results.levels;
    vals{i} = F2(levels, :);
    
    % calculate x and z:
    % convert B from MHz to Gauss
    B = B / p.Results.gF; % gF uB in MHz/G
    x = ((1-cos(theta)) * (B / qdrpGrad).^2).^0.5 * 1e4;
    z = (cos(theta) * (B / qdrpGrad).^2 / 4).^0.5 * 1e4;
    
    % upper hemisphere (to check properly)
    if ~isreal(z) 
        z = -abs(z);
        x = ((1+cos(theta)) * (B / qdrpGrad).^2).^0.5 * 1e4;
    end
    
    xs{i} = x;
    zs{i} = z;
    
    if p.Results.showGraphs
        hold on; plot3(x,-z,F2(levels, :), 'k.'); hold on; pause(0.01);
        view([-45 50]); %axis equal
        xlabel('X (\mum)');
        ylabel('Z (\mum)');
        title('Shell trap');
    end
    
    fprintf('Calculating theta=%2.2f...\n', theta)
end

fprintf('Evaluation complete.\n')

% Save the generated data to a binary file. Data is saved in strip line
% blocks. Each stripline consists of the following format:
% THETA NUMBER_POINT X Z V

    function saveStrip(fH, theta, N, x, z, v)
        fwrite(fH, theta, 'double');
        fwrite(fH, N, 'double');
        fwrite(fH, x, 'double');
        fwrite(fH, z, 'double');
        
        % write each level together.
        fwrite(fH, size(v, 1), 'uint32');
        fwrite(fH, v', 'double');
    end

fH = fopen(p.Results.filename, 'w');
oc = onCleanup(@() fclose(fH));

fwrite(fH, size(RF, 1), 'uint32');
fwrite(fH, RF, 'double');
fwrite(fH, size(Rabi, 1), 'uint32');
fwrite(fH, Rabi, 'double');
fwrite(fH, qdrpGrad, 'double');

for i=1:length(thetas)
    saveStrip(fH, thetas(i), length(xs{i}), xs{i}, zs{i}, vals{i});
end

end

