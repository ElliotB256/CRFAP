% Try to generate spectrum of noise...

RF = 1;
gRF = 0.3;

rates = [];
RFps=0:0.05:2;
for RFp=RFps
    rates(end+1) = getRate(RF, RFp, gRF, 0.1);
    fprintf('RFp=%.3f\n', RFp)
end
plot(RFps,rates);