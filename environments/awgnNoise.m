function x = awgnNoise(xClean, noiseCfg)
%AWGNNOISE Add additive white Gaussian noise at the requested SNR.

signalPower = mean(xClean.^2);
noisePower = signalPower / 10^(noiseCfg.SNR_dB/10);
noise = sqrt(noisePower) * randn(size(xClean));

x = xClean + noise;
end
