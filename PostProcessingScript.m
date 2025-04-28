%% Post Processing 
clear all;
close all;
clc;

%% Load Data and Align Samples

% Radar Parameters
Fs = 60e6;
Ts = 1/Fs;
BW = 30e6;
frequencyStart = 2.45e9 - BW/2;
frequencyEnd = 2.45e9 + BW/2;
fo = (frequencyStart + frequencyEnd)/2;

tau =10e-6;
u = BW / tau;

sti = load_sc16q11_MIMO('R:\Temp\transmit.sc16q11', 2);
sri = load_sc16q11_MIMO('R:\Temp\receive.sc16q11', 2);

% Align samples 
if length(sti) > length(sri)
    sti = sti(1:size(sri,1), :); % Trim sti to match sri row count
    t = (0:length(sri)-1) / Fs; % Time vector
else
    sri = sri(1:size(sti,1), :); % Trim sri to match sti row count
    t = (0:length(sti)-1) / Fs; % Time vector
end

%% Plot TX vs RX 
% Real and Imaginary Parts of Chirp vs. Received Signal
figure('Color', [1 1 1]);
subplot(2,2,1);
plot(t, real(sti(:,1)), 'b', t, imag(sti(:,1)), 'r');
title('Chirp Signal - TX1');
xlabel('Time (s)');
xlim([0 10e-6]);
ylabel('Amplitude');
ylim([-1 1]);
legend('Real Part', 'Imaginary Part');
grid on;

subplot(2,2,2);
plot(t, real(sti(:,2)), 'b', t, imag(sti(:,2)), 'r');
title('Chirp Signal - TX2');
xlabel('Time (s)');
xlim([0 10e-6]);
ylabel('Amplitude');
ylim([-1 1]);
legend('Real Part', 'Imaginary Part');
grid on;

subplot(2,2,3);
plot(t, real(sri(:,1)), 'b', t, imag(sri(:,1)), 'r');
title('Chirp Signal - RX1');
xlabel('Time (s)');
xlim([0 10e-6]);
ylabel('Amplitude');
ylim([-1 1]);
legend('Real Part', 'Imaginary Part');
grid on;

subplot(2,2,4);
plot(t, real(sri(:,2)), 'b', t, imag(sri(:,2)), 'r');
title('Chirp Signal - RX2');
xlabel('Time (s)');
xlim([0 10e-6]);
ylabel('Amplitude');
ylim([-1 1]);
legend('Real Part', 'Imaginary Part');
grid on;

%% FFT TX vs RX
% Chirp vs. Received 

N = length(sti(:,1));                   % Number of samples
X = fftshift(abs(fft(sti(:,1), N)));    % Compute the shifted FFT
f = linspace(-Fs/2,Fs/2,N);             % Frequency axis centered -Fs/2 to Fs/2
f = f + (frequencyStart);               % Shift to carrier frequency (Hz)

figure('Color', [1 1 1]);
subplot(2,2,1);
plot(f / 1e9, abs(X) / max(abs(X))); % Normalize and scale to GHz
title('FFT Spectrum - TX1');
xlabel('Frequency (GHz)');
xlim([2.435 2.465]);
ylabel('Normalized Magnitude');
ylim([0 1]);
grid on;

N = length(sti(:,2));                   % Number of samples
X = fftshift(abs(fft(sti(:,2), N)));    % Compute the shifted FFT
f = linspace(-Fs/2,Fs/2,N);             % Frequency axis centered -Fs/2 to Fs/2
f = f + (frequencyStart);               % Shift to carrier frequency (Hz)

subplot(2,2,3);
plot(f / 1e9, abs(X) / max(abs(X))); % Normalize and scale to GHz
title('FFT Spectrum - RX1');
xlabel('Frequency (GHz)');
xlim([2.435 2.465]);
ylabel('Normalized Magnitude');
ylim([0 1]);
grid on;

N = length(sri(:,1));                   % Number of samples
X = fftshift(abs(fft(sri(:,1), N)));    % Compute the shifted FFT
f = linspace(-Fs/2,Fs/2,N);             % Frequency axis centered -Fs/2 to Fs/2
f = f + (frequencyStart);               % Shift to carrier frequency (Hz)

subplot(2,2,2);
plot(f / 1e9, abs(X) / max(abs(X))); % Normalize and scale to GHz
title('FFT Spectrum - TX2');
xlabel('Frequency (GHz)');
xlim([2.435 2.465]);
ylabel('Normalized Magnitude');
ylim([0 1]);
grid on;

N = length(sri(:,2));                   % Number of samples
X = fftshift(abs(fft(sri(:,2), N)));    % Compute the shifted FFT
f = linspace(-Fs/2,Fs/2,N);             % Frequency axis centered -Fs/2 to Fs/2
f = f + (frequencyStart);               % Shift to carrier frequency (Hz)

subplot(2,2,4);
plot(f / 1e9, abs(X) / max(abs(X))); % Normalize and scale to GHz
title('FFT Spectrum - RX2');
xlabel('Frequency (GHz)');
xlim([2.435 2.465]);
ylabel('Normalized Magnitude');
ylim([0 1]);
grid on;

%% SYSTEM DELAY Estimation 
% Extract real parts of TX and RX channels
txCH1 = real(sti(:, 1));
txCH2 = real(sti(:, 2));
rxCH1 = real(sri(:, 1));
rxCH2 = real(sri(:, 2));

%%%%%%%%% COARSE SYSTEM DELAY %%%%%%%%%%%%%%%%%%%%%%%%%% 
% Cross-Correlation to Determine Start Indices
[corrCH1, lagsCH1] = xcorr(rxCH1, txCH1);
[~, maxIdxCH1] = max(abs(corrCH1));
startIdxCH1_Coarse = lagsCH1(maxIdxCH1);

[corrCH2, lagsCH2] = xcorr(rxCH2, txCH2);
[~, maxIdxCH2] = max(abs(corrCH2));
startIdxCH2_Coarse = lagsCH2(maxIdxCH2);

%%%%%%%%% FINE SYSTEM DELAY (Using Interpolation) %%%%%%
% Apply Low-Pass Filtering Before Cross-Correlation
lp_cutoff = BW/2; % Cutoff frequency = half of signal bandwidth
sti_filtered = lowpass(real(sti), lp_cutoff, Fs) + 1i * lowpass(imag(sti), lp_cutoff, Fs);
sri_filtered = lowpass(real(sri), lp_cutoff, Fs) + 1i * lowpass(imag(sri), lp_cutoff, Fs);

interp_factor = 300; 
t_original = 1:length(sti_filtered(:,1));
t_interpolated = linspace(1, length(sti_filtered(:,1)), interp_factor * length(sti_filtered(:,1)));
sti_interp = interp1(t_original, real(sti_filtered), t_interpolated, 'spline', 'extrap');
sri_interp = interp1(t_original, real(sri_filtered), t_interpolated, 'spline', 'extrap');

% Cross-correlation of up-sampled signals
[corrCH1, lagsCH1] = xcorr(sri_interp(:,1), sti_interp(:,1));
[~, maxidxCH1] = max(abs(corrCH1));
startIdxCH1_Fine = lagsCH1(maxidxCH1)/interp_factor; % Sub-Sample Delay Estimation

[corrCH2, lagsCH2] = xcorr(sri_interp(:,2), sti_interp(:,2));
[~, maxidxCH2] = max(abs(corrCH2));
startIdxCH2_Fine = lagsCH2(maxidxCH2)/interp_factor; % Sub-Sample Delay Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Known System Delay Compensation (in samples)
systemDelaySamplesCH1 = 42;
systemDelaySamplesCH2 = 41;
systemDelayCH1_Coarse = systemDelaySamplesCH1 / Fs;
systemDelayCH2_Coarse = systemDelaySamplesCH2 / Fs;

% Coarse Delay Estimation
tdCH1_Coarse = startIdxCH1_Coarse / Fs;
tdCH2_Coarse = startIdxCH2_Coarse / Fs;
tdCH1_Fine = startIdxCH1_Fine / (Fs * interp_factor);
tdCH2_Fine = startIdxCH2_Fine / (Fs * interp_factor);

% Phase Delay Estimation
phase_tx1 = angle(sti(:,1));
phase_tx2 = angle(sti(:,2));
phase_rx1 = angle(sri(:,1));
phase_rx2 = angle(sri(:,2));
phi_diff = unwrap((phase_rx1 - phase_tx1) + (phase_rx2 - phase_tx2)) / 2;
phi_diff_deg = rad2deg(mean(phi_diff));
td_Phase = mean(phi_diff) / (2 * pi * fo);
fprintf('\nCoarse Delay Estimation Results (Hardcoded System Delay):\n');
fprintf('----------------------------------------------------------\n');
fprintf('StartIdx CH1 (Coarse): %.3f samples | Delay: %.3f ns\n', startIdxCH1_Coarse, tdCH1_Coarse * 1e9);
fprintf('StartIdx CH2 (Coarse): %.3f samples | Delay: %.3f ns\n', startIdxCH2_Coarse, tdCH2_Coarse * 1e9);
fprintf('StartIdx CH1   (Fine): %.3f samples | Delay: %.3f ns\n', startIdxCH1_Fine, tdCH1_Fine * 1e9);
fprintf('StartIdx CH2   (Fine): %.3f samples | Delay: %.3f ns\n', startIdxCH2_Fine, tdCH2_Fine * 1e9);
fprintf('Phase Difference:      %.2f deg   | Time Delay: %.3f ns\n', phi_diff_deg, td_Phase * 1e9);

%% Final Delay Estimation (Subtract System Delay AND Phase Correction)
tdCH1_System = tdCH1_Coarse - systemDelayCH1_Coarse - td_Phase;
tdCH2_System = tdCH2_Coarse - systemDelayCH2_Coarse - td_Phase;
tdCoarse_avg = (tdCH1_System + tdCH2_System)/2;
tdFine_avg = (tdCH1_Fine + tdCH2_Fine)/2;
tdCH1_total = tdCH1_Coarse + tdCH1_Fine;
tdCH2_total = tdCH2_Coarse + tdCH2_Fine;
%td = (tdCoarse_avg + tdFine_avg)/2; % This line is for calibrating with 1ft sample
td = tdCoarse_avg + tdFine_avg;

% Range Calculation
c = physconst('LightSpeed');
rangeCH1_Coarse = tdCH1_System * c / 2;
rangeCH2_Coarse = tdCH2_System * c / 2;
range_Coarse = tdCoarse_avg * c / 2;
rangeCoarse_ft = range_Coarse * 3.28084;
rangeCH1_Fine = tdCH1_Fine * c / 2;
rangeCH2_Fine = tdCH2_Fine * c / 2;
range_Fine = tdFine_avg * c / 2;
rangeFine_ft = range_Fine * 3.28084;
range_td = td * c / 2;
range_td_ft = range_td * 3.28084;

% Display Results
fprintf('\nChannel 1 - Range and Delay:\n');
fprintf('----------------------------------------------------------\n');
fprintf('CH1 Delay (Corrected): %.3f ns | Range: %.3f m (%.3f ft)\n', tdCH1_Coarse*1e9, rangeCH1_Coarse, rangeCH1_Coarse*3.28084);
fprintf('CH2 Delay (Corrected): %.3f ns | Range: %.3f m (%.3f ft)\n', tdCH2_Coarse*1e9, rangeCH2_Coarse, rangeCH2_Coarse*3.28084);
fprintf('Averaged Delay: %.3f ns | Estimated Range: %.3f m (%.3f ft)\n', tdCoarse_avg*1e9, range_Coarse, rangeCoarse_ft);
fprintf('----------------------------------------------------------\n');
fprintf('\nChannel 2 - Range and Delay:\n');
fprintf('----------------------------------------------------------\n');
fprintf('CH1 Delay (Corrected): %.3f ns | Range: %.3f m (%.3f ft)\n', tdCH1_Fine*1e9, rangeCH1_Fine, rangeCH1_Fine*3.28084);
fprintf('CH2 Delay (Corrected): %.3f ns | Range: %.3f m (%.3f ft)\n', tdCH2_Fine*1e9, rangeCH2_Fine, rangeCH2_Fine*3.28084);
fprintf('Averaged Delay: %.3f ns | Estimated Range: %.3f m (%.3f ft)\n', tdFine_avg*1e9, range_Fine, rangeFine_ft);
fprintf('----------------------------------------------------------\n');

fprintf('\nFinal Results - Range and Delay:\n');
fprintf('----------------------------------------------------------\n');
fprintf('Final Delay: %.3f ns | Estimated Range: %.3f m (%.3f ft)\n', td*1e9, range_td, range_td_ft);

%% Obtain Beat Frequencies
% The beat signal is extracted from the complex signals. By taking the
% conjugate of the received signal the frequency difference can be
% extracted.
beatSignal1 = sti(:,1) .* conj(sri(:,1)); 
beatSignal2 = sti(:,2) .* conj(sri(:,2));

%td_residual = Ts - abs(td_Phase) - tdCH1_Fine;        % Attemots to account for residual fine delay errors
%tau_system1 = systemDelaySamplesCH1/Fs - td_residual;
%tau_system2 = systemDelaySamplesCH2/Fs - td_residual;

tau_system1 = systemDelaySamplesCH1/Fs - 1.1295e-08; % Most accurate 
tau_system2 = systemDelaySamplesCH2/Fs - 1.1295e-08;

%tau_system1 = (startIdxCH1_Coarse) / Fs + tdCH1_Fine; % These two line are needed when calibrating with 1 ft loopback test 
%tau_system2 = (startIdxCH2_Coarse) / Fs + tdCH2_Fine; %

t = (0:length(beatSignal1)-1).' / Fs;

% Time-domain delay removal
beatSignal1 = beatSignal1 .* exp(-1j * 2 * pi * u * t * tau_system1);
beatSignal2 = beatSignal2 .* exp(-1j * 2 * pi * u * t * tau_system2);



%% FFT Processing
% The Fast Fourier Transform is performed on the downsmapled beat signal. 
% The result is the beat frequency fb which is used later in the script to obtain the range and velocity.

% FFT
Nfft = 2^24; % Resolution (24 is rather high but acheives very good results)
f_beat = (-Nfft/2:Nfft/2-1) * (Fs / Nfft); % Frequency vector for plotting
beatFFT1 = fftshift(fft(beatSignal1, Nfft));
beatFFT2 = fftshift(fft(beatSignal2, Nfft));

% Plot Spectrums
figure('Color', [1 1 1]);
subplot(2,1,1);
plot(f_beat/1e6, 20*log10(abs(beatFFT1)));                   % Plot in MHz
title('Beat Frequency Spectrum - Element 1');
xlabel('Frequency (MHz)');
xlim([-5 5]);
ylabel('Magnitude (dB)');
ylim([-10 60]);
grid on;
subplot(2,1,2);
plot(f_beat/1e6, 20*log10(abs(beatFFT2)));                  
title('Beat Frequency Spectrum - Element 2');
xlabel('Frequency (MHz)');
xlim([-5 5]);
ylabel('Magnitude (dB)');
ylim([-10 60]);
grid on;

%% MUSIC Algorithm for DOA Estimation
% based on B. Allen and M. Ghavami, Adaptive Array Systems. John Wiley & Sons, 2006.
% One can estimate the steering vectors associated with the received signals by finding 
% the steering vectors which are most nearly orthogonal to the eigenvectors associated with the eigenvalues of R. 

% Antenna Array parameters
c = physconst('LightSpeed'); % Speed of light
lambda = c/fo; % wavelength
d = lambda/2;  % Half-wavelength spacing
num_elements = 2; % Two elements in TX and RX array
element_positions = (0:num_elements-1) * d; % positions of antenna using half-wavelength spacing

% Set length
snapshotLength = 50;                               % Set snapshot length to 50

% Set start and end index values
if startIdxCH1_Coarse > startIdxCH2_Coarse
startIdx = startIdxCH1_Coarse; % Set start index to 1
else 
    startIdx = startIdxCH2_Coarse;
end 
endIdx = startIdx + snapshotLength - 1;                 % Calculate end index

% Create array data matrix using complex samples
X = zeros(2, snapshotLength);                           % Initialize a 2 row matrix for sample data
X(1,:) = beatSignal1(startIdx:endIdx);       % Store downsampled first beat signal in first row of X
X(2,:) = beatSignal2(startIdx:endIdx);       % Store downsampled second beat signal in second row of X

R = X * X' / snapshotLength;                            % Calculate spatial correlation matrix

% Eigendecomposition
[V, D] = eig(R);                                        % Perform eigendecomposition on matrix R
[D, idx] = sort(diag(D), 'descend');                    % Sort eigenvalues in descending order
V = V(:, idx);                                          % Reorder eigenvectors to match sorted eigenvalues

% Separate signal and noise subspaces
targets = 1;                                            % Single target
En = V(:, targets+1:end);                               % Noise subspace

% Scan angles and compute MUSIC spectrum
thetaScan = -90:0.1:90;                                 % Define scan angles
Pmusic = zeros(size(thetaScan));                        % Initialize MUSIC spectrum array

for i = 1:length(thetaScan)
    
    % Steering vector
    a = exp(-1j*2*pi*d*sind(thetaScan(i))/lambda * (0:num_elements-1).');  
   
    % MUSIC spectrum
    den = a' * (En * En') * a;     

    % Calculate MUSIC spectrum
    Pmusic(i) = 1 / abs(den);   

end

% Normalize spectrum
Pmusic = Pmusic / max(Pmusic);

% Plot MUSIC spectrum in dB's
figure('Color', [1 1 1]);
plot(thetaScan, 10*log10(Pmusic)); 
title('MUSIC Spectrum for DOA Estimation');
xlabel('Angle (degrees)');
xlim([-90 90]);
ylabel('Normalized Power (dB)');
ylim([-80 20]);
grid on;

% Find peaks in MUSIC spectrum
[pks, locs] = findpeaks(10*log10(Pmusic), thetaScan);   % Find peaks in the MUSIC spectrum

if ~isempty(locs)
    AoA = locs(1);           % Take first peak
    pk_val = pks(1);

    % Add peak annotation to the plot
    hold on;
    plot(AoA, pk_val, 'r.', 'MarkerSize', 20);
    text(AoA, pk_val, [' ' num2str(AoA, '%.1f') '°'], ...
        'VerticalAlignment', 'bottom', 'Color', 'r');
    hold off;
else
    % No directionality present
    AoA = NaN;
    pk_val = NaN;
    warning('No peaks found in MUSIC spectrum. Likely due to lack of spatial diversity in signal.');
end

%% Range and Velocity Calculations

% Split FFT into positive and negative halves
beatFFT1_pos = abs(beatFFT1(Nfft/2+1:end));  % positive freqs
beatFFT1_neg = abs(beatFFT1(1:Nfft/2));      % negative freqs
[max_val_pos1, idxPos1] = max(beatFFT1_pos);
[max_val_neg1, idxNeg1] = max(beatFFT1_neg);
beatFFT2 = fftshift(fft(beatSignal2, Nfft));
beatFFT2_pos = abs(beatFFT2(Nfft/2+1:end));
beatFFT2_neg = abs(beatFFT2(1:Nfft/2));
[max_val_pos2, idxPos2] = max(beatFFT2_pos);
[max_val_neg2, idxNeg2] = max(beatFFT2_neg);

% Map indices to frequency bins
BeatFrequencyUpper1 = f_beat(Nfft/2 + idxPos1);
BeatFrequencyLower1 = f_beat(idxNeg1);
BeatFrequencyUpper2 = f_beat(Nfft/2 + idxPos2);
BeatFrequencyLower2 = f_beat(idxNeg2);

if max_val_pos1 >= max_val_neg1
    fb1 = abs(f_beat(Nfft/2 + idxPos1));
else
    fb1 = abs(f_beat(idxNeg1));
end

if max_val_pos2 >= max_val_neg2
    fb2 = abs(f_beat(Nfft/2 + idxPos2));
else
    fb2 = abs(f_beat(idxNeg2));
end

%fb = fb1;
fb = (fb1 + fb2) / 2;

fprintf('Upper Beat Frequency: %.2f Hz\n', BeatFrequencyUpper1);
fprintf('Lower Beat Frequency: %.2f Hz\n', BeatFrequencyLower1);
fprintf('Upper Beat Frequency: %.2f Hz\n', BeatFrequencyUpper2);
fprintf('Lower Beat Frequency: %.2f Hz\n', BeatFrequencyLower2);

tau_from_fb = fb / u;              % time delay (seconds)
range_fb = (c * tau_from_fb) / 2;  % one-way range (meters)
range_fb_ft = range_fb * 3.28084; % convert to feet

range = (range_fb + range_td)/2;
range_ft = range * 3.28084; % convert to feet

fprintf('\nRange from Beat Frequency (Delay from Frequency Domain):\n');
fprintf('----------------------------------------------------------\n');
fprintf('Beat Frequency: %.3f kHz\n', fb * 1e-3);
fprintf('Estimated Time Delay: %.3f ns\n', tau_from_fb * 1e9);
fprintf('Estimated Range: %.3f m (%.3f ft)\n', range_fb, range_fb_ft);

%% Doppler and Velocity Estimation

% Choose the stronger peak
if fb == BeatFrequencyLower1
    fb_doppler = abs(BeatFrequencyUpper1);
else
    fb_doppler = abs(BeatFrequencyLower1);
end

%fb_doppler = abs(BeatFrequencyLower);  % Use upper for Doppler sense
fd = fb_doppler * 2 * pi;
velocity = (fd * lambda) / (2 * cosd(AoA));

fprintf('\nRadar Measurement Results:\n');
fprintf('----------------------------\n');
fprintf('Estimated Beat Frequency: %.2f kHz\n', fb * 1e-3);
fprintf('Estimated Doppler Frequency: %.2f Hz\n', fd);
fprintf('Estimated Range: %.2f m (%.2f ft)\n', range, range_ft);
fprintf('Estimated Velocity: %.2f m/s\n', velocity);
fprintf('Estimated Angle: %.2f degrees\n', AoA);
fprintf('----------------------------\n');

td_trusted = td;  % Your delay estimate after correction
fb_from_td = u * td_trusted;
fprintf('Expected Beat Frequency from Delay: %.2f kHz\n', fb_from_td * 1e-3);

error_kHz = abs(fb - fb_from_td) / 1e3;
fprintf('Beat frequency error vs. delay-based: %.3f kHz\n', error_kHz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% AoA = detected angle (degrees)
% range_ft = detected range (feet)
% thetaScan = -90:0.5:90;

% Build a tight range grid around the target
range_axis = linspace(0, 40, 1000); % 1000 spaced between 0 and 40 feet
[ThetaGrid, RangeGrid] = meshgrid(thetaScan, range_axis); % build a meshgrid using thetascan and range_axis
TargetMap = zeros(size(RangeGrid)); % stores intensity values 

kernelSize = 51;           % Will be used to fill in a heatmap circle. Odd values will allow for an accurate center. 
halfSize = floor(kernelSize / 2); % Cutting the 
sigma = kernelSize / 9;   % controls the spread of intensity
[x, y] = meshgrid(-halfSize:halfSize, -halfSize:halfSize); 
intensityKernel = exp(-(x.^2 + y.^2) / (2 * sigma^2)); % Gaussian kernel function
intensityKernel = intensityKernel / max(intensityKernel(:));  % normalize peak to 1

% Find location for time-domain estimated range
[~, angleIdx1] = min(abs(thetaScan - AoA)); % Location on x axis
[~, rangeIdx1] = min(abs(range_axis - range_td_ft)); % Location on y axis

% Find location for frequency-domain estimated range
[~, rangeIdx2] = min(abs(range_axis - range_fb_ft));
[~, angleIdx2] = min(abs(thetaScan - AoA)); % Same AoA

% Average estimated range
[~, rangeIdx3] = min(abs(range_axis - range_ft));
[~, angleIdx3] = min(abs(thetaScan - AoA)); % Same AoA

% Place both kernels
for i = -halfSize:halfSize
    for j = -halfSize:halfSize
        
        % Time-domain kernel
        ri1 = rangeIdx1 + i;
        ai1 = angleIdx1 + j;
        if ri1 > 0 && ri1 <= size(TargetMap, 1) && ai1 > 0 && ai1 <= size(TargetMap, 2)
            TargetMap(ri1, ai1) = TargetMap(ri1, ai1) + intensityKernel(i + halfSize + 1, j + halfSize + 1); % weight = 0.8
        end

        % Frequency-domain kernel
        ri2 = rangeIdx2 + i;
        ai2 = angleIdx2 + j;
        if ri2 > 0 && ri2 <= size(TargetMap, 1) && ai2 > 0 && ai2 <= size(TargetMap, 2)
            TargetMap(ri2, ai2) = TargetMap(ri2, ai2) + intensityKernel(i + halfSize + 1, j + halfSize + 1); % weight = 0.6
        end

        % Frequency-domain kernel
        ri3 = rangeIdx3 + i;
        ai3 = angleIdx3 + j;
        if ri3 > 0 && ri3 <= size(TargetMap, 1) && ai3 > 0 && ai3 <= size(TargetMap, 2)
            TargetMap(ri3, ai3) = TargetMap(ri3, ai3) + intensityKernel(i + halfSize + 1, j + halfSize + 1); % weight = 0.6
        end

    end
end

% Plot
figure('Color', [1 1 1]);
imagesc(thetaScan, range_axis, TargetMap); % Plot as a color image
xlabel('Angle of Arrival (degrees)');
ylabel('Range (ft)');
title('Range vs AoA Heatmap');
axis xy; % Y axis goes from bottom to top
colormap(jet); % color - blue(low) to red(high)
colorbar; % shows intensity scale
clim([0 1]); % Forces max intensity
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
