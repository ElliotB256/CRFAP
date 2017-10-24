function Sample(instance)
%SAMPLE Calculate eigenenergies at each B.

% First, create a number of line sampler objects. Each is configured
% according to properties of the instance.

    function lineSampler = createSubsampler(instance)
       lineSampler = AP.Sampler.LineSampler(instance.APCalculator);
       lineSampler.Gamma = instance.Gamma;
       lineSampler.StartB = instance.StartB;
       lineSampler.MeshIterations = instance.MeshIterations;
       lineSampler.QuadGrad = instance.QuadGrad;
       lineSampler.Sort = instance.SortLines;
       lineSampler.Verbose = instance.Verbose;
    end

thetas = linspace(instance.ThetaRange(1), instance.ThetaRange(2), instance.RayNumber);
instance.LineSamplers = cell(1, length(thetas));
for i=1:length(thetas)
    ls = createSubsampler(instance);
    ls.Theta = thetas(i);
    instance.LineSamplers{i} = ls;
end

% Each line sampler is now sampled in turn.
fprintf('Created %d total line samplers.\n', length(instance.LineSamplers))
for i=1:length(instance.LineSamplers)
   s = instance.LineSamplers{i};
   fprintf('Sampling Line %d of %d:\n', i, length(instance.LineSamplers))
   s.Sample();
end

% Set instance to clean
instance.Dirty = 0;

end