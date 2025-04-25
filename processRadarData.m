function [td, range_td_ft, tau_from_fb, range_fb_ft, range_ft, ...
          startIdxCH1_Coarse, startIdxCH2_Coarse] = ...
          processRadarData(txFile, rxFile, systemDelaySamplesCH1, systemDelaySamplesCH2)

% Constants
Fs = 60e6;
BW = 30e6;
tau = 10e-6;
u = BW / tau;
c = physconst('LightSpeed');
fo = 2.45e9;

% Load signals
sti = load_sc16q11_MIMO(txFile, 2);
sri = load_sc16q11_MIMO(rxFile, 2);
minLen = min(size(sti,1), size(sri,1));
sti = sti(1:minLen,:);
sri = sri(1:minLen,:);
t = (0:minLen-1).'/Fs;

% Coarse delay via cross-correlation
[cCH1, lCH1] = xcorr(real(sri(:,1)), real(sti(:,1)));
[~, idxCH1] = max(abs(cCH1));
startIdxCH1_Coarse = lCH1(idxCH1);

[cCH2, lCH2] = xcorr(real(sri(:,2)), real(sti(:,2)));
[~, idxCH2] = max(abs(cCH2));
startIdxCH2_Coarse = lCH2(idxCH2);

% Fine delay via interpolation
sti_filt = lowpass(real(sti), BW/2, Fs) + 1j * lowpass(imag(sti), BW/2, Fs);
sri_filt = lowpass(real(sri), BW/2, Fs) + 1j * lowpass(imag(sri), BW/2, Fs);
interp_factor = 300;
t_orig = 1:minLen;
t_interp = linspace(1, minLen, interp_factor * minLen);
sti_interp = interp1(t_orig, real(sti_filt), t_interp, 'spline', 'extrap');
sri_interp = interp1(t_orig, real(sri_filt), t_interp, 'spline', 'extrap');

[cCH1_f, lCH1_f] = xcorr(sri_interp(:,1), sti_interp(:,1));
[~, idxCH1_f] = max(abs(cCH1_f));
startIdxCH1_Fine = lCH1_f(idxCH1_f) / interp_factor;

[cCH2_f, lCH2_f] = xcorr(sri_interp(:,2), sti_interp(:,2));
[~, idxCH2_f] = max(abs(cCH2_f));
startIdxCH2_Fine = lCH2_f(idxCH2_f) / interp_factor;

% Phase delay
phi_diff = unwrap(angle(sri(:,1)) - angle(sti(:,1)) + angle(sri(:,2)) - angle(sti(:,2))) / 2;
td_Phase = mean(phi_diff) / (2*pi*fo);

% Delay (time domain)
tdCH1_Coarse = startIdxCH1_Coarse / Fs;
tdCH2_Coarse = startIdxCH2_Coarse / Fs;
tdCH1_Fine = startIdxCH1_Fine / (Fs * interp_factor);
tdCH2_Fine = startIdxCH2_Fine / (Fs * interp_factor);

sysCH1 = systemDelaySamplesCH1 / Fs;
sysCH2 = systemDelaySamplesCH2 / Fs;

tdCH1 = tdCH1_Coarse - sysCH1 - td_Phase;
tdCH2 = tdCH2_Coarse - sysCH2 - td_Phase;

tdCoarse = (tdCH1 + tdCH2) / 2;
tdFine = (tdCH1_Fine + tdCH2_Fine) / 2;
td = tdCoarse + tdFine;

range_td_ft = td * c / 2 * 3.28084;

% Frequency domain estimation
beat1 = sti(:,1) .* conj(sri(:,1));
beat2 = sti(:,2) .* conj(sri(:,2));
tau1 = systemDelaySamplesCH1/Fs - 1.1295e-08;
tau2 = systemDelaySamplesCH2/Fs - 1.1295e-08;
beat1 = beat1 .* exp(-1j * 2 * pi * u * t * tau1);
beat2 = beat2 .* exp(-1j * 2 * pi * u * t * tau2);

N = 2^24;
f = (-N/2:N/2-1) * Fs / N;
B1_fft = fftshift(fft(beat1, N));
B1_pos = abs(B1_fft(N/2+1:end));
B1_neg = abs(B1_fft(1:N/2));
[~, i_pos] = max(B1_pos);
[~, i_neg] = max(B1_neg);

if B1_pos(i_pos) >= B1_neg(i_neg)
    fb = abs(f(N/2 + i_pos));
else
    fb = abs(f(i_neg));
end

tau_from_fb = fb / u;
range_fb_ft = (c * tau_from_fb / 2) * 3.28084;
range_ft = ((range_td_ft / 3.28084 + range_fb_ft / 3.28084)/2) * 3.28084;

end
