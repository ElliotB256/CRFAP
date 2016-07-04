function b = Unique(a, thresh)
% UNIQUE Gets unique elements in a within threshold (works for floats)

b = [];
for i=1:length(a)
    if any(abs(b-a(i)) < thresh)
        continue;
    else
        b(end+1) = a(i);
end


end