% Test we are generating uthe right results from the dressed state basis
Bs = 0.0:0.05:2;

figure(1);

L = [];

for B = Bs
    [E,Ev,~] = dressedStateBasis(B, [2 1.5], [6 1], [0.2 0]);
    L(end+1, :) = diag(Ev);
    hold on;
    plot(B*ones(size(Ev)), Ev, '.k');
    hold off;
end

plot(Bs,sort(L,2));