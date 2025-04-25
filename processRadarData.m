function [td, range_td_ft, tau_from_fb, range_fb_ft, range_ft, ...
          startIdxCH1_Coarse, startIdxCH2_Coarse, ...
          systemDelaySamplesCH1, systemDelaySamplesCH2, ...
          thetaScan, range_axis, TargetMap, ...
          f_rx, rxFFT1, rxFFT2, ...
          sti, sri, Fs, frequencyStart, frequencyEnd, fo, t, ...
                  beatFFT1, beatFFT2, f_beat, AoA] = ...
          processRadarData(txFile, rxFile, systemDelaySamplesCH1, systemDelaySamplesCH2)

% Radar Parameters
Fs = 60e6;
BW = 30e6;
frequencyStart = 2.45e9 - BW/2;
frequencyEnd = 2.45e9 + BW/2;
fo = (frequencyStart + frequencyEnd)/2;
tau = 10e-6;
u = BW / tau;
c = physconst('LightSpeed');

% Load Signals
sti = load_sc16q11_MIMO(txFile, 2);
sri = load_sc16q11_MIMO(rxFile, 2);

% Align sample lengths
t = (0:min(length(sti), length(sri))-1) / Fs;
sti = sti(1:length(t), :);
sri = sri(1:length(t), :);

% Coarse Delay
[corrCH1, lagsCH1] = xcorr(real(sri(:,1)), real(sti(:,1)));
[~, maxIdxCH1] = max(abs(corrCH1));
startIdxCH1_Coarse = lagsCH1(maxIdxCH1);

[corrCH2, lagsCH2] = xcorr(real(sri(:,2)), real(sti(:,2)));
[~, maxIdxCH2] = max(abs(corrCH2));
startIdxCH2_Coarse = lagsCH2(maxIdxCH2);

% Fine Delay (interpolated)
lp_cutoff = BW/2;
sti_filtered = lowpass(real(sti), lp_cutoff, Fs) + 1i * lowpass(imag(sti), lp_cutoff, Fs);
sri_filtered = lowpass(real(sri), lp_cutoff, Fs) + 1i * lowpass(imag(sri), lp_cutoff, Fs);
interp_factor = 300;
t_original = 1:length(sti_filtered);
t_interpolated = linspace(1, length(sti_filtered), interp_factor * length(sti_filtered));
sti_interp = interp1(t_original, real(sti_filtered), t_interpolated, 'spline', 'extrap');
sri_interp = interp1(t_original, real(sri_filtered), t_interpolated, 'spline', 'extrap');

[corrCH1, lagsCH1] = xcorr(sri_interp(:,1), sti_interp(:,1));
[~, maxidxCH1] = max(abs(corrCH1));
startIdxCH1_Fine = lagsCH1(maxidxCH1)/interp_factor;

[corrCH2, lagsCH2] = xcorr(sri_interp(:,2), sti_interp(:,2));
[~, maxidxCH2] = max(abs(corrCH2));
startIdxCH2_Fine = lagsCH2(maxidxCH2)/interp_factor;

% Phase delay
phi_diff = unwrap((angle(sri(:,1)) - angle(sti(:,1))) + (angle(sri(:,2)) - angle(sti(:,2)))) / 2;
td_Phase = mean(phi_diff) / (2 * pi * fo);

% Delay Calculations
systemDelayCH1_Coarse = systemDelaySamplesCH1 / Fs;
systemDelayCH2_Coarse = systemDelaySamplesCH2 / Fs;
tdCH1_Coarse = startIdxCH1_Coarse / Fs;
tdCH2_Coarse = startIdxCH2_Coarse / Fs;
tdCH1_Fine = startIdxCH1_Fine / (Fs * interp_factor);
tdCH2_Fine = startIdxCH2_Fine / (Fs * interp_factor);

tdCH1_System = tdCH1_Coarse - systemDelayCH1_Coarse - td_Phase;
tdCH2_System = tdCH2_Coarse - systemDelayCH2_Coarse - td_Phase;
tdCoarse_avg = (tdCH1_System + tdCH2_System)/2;
tdFine_avg = (tdCH1_Fine + tdCH2_Fine)/2;
%td = (tdCoarse_avg + tdFine_avg)/2; % This line is for calibrating with 1ft sample
td = tdCoarse_avg + tdFine_avg;

range_td = td * c / 2;
range_td_ft = range_td * 3.28084;

% Beat Frequency
beatSignal1 = sti(:,1) .* conj(sri(:,1));
beatSignal2 = sti(:,2) .* conj(sri(:,2));

%td_residual = Ts - abs(td_Phase) - tdCH1_Fine;        % Attemots to account for residual fine delay errors
%tau_system1 = systemDelaySamplesCH1/Fs - td_residual;
%tau_system2 = systemDelaySamplesCH2/Fs - td_residual;

tau_system1 = systemDelaySamplesCH1/Fs - 1.1295e-08; % Most accurate 
tau_system2 = systemDelaySamplesCH2/Fs - 1.1295e-08;

%tau_system1 = (startIdxCH1_Coarse) / Fs + tdCH1_Fine; % These two line are needed when calibrating with 1 ft loopback test 
%tau_system2 = (startIdxCH2_Coarse) / Fs 

t = (0:length(beatSignal1)-1).' / Fs;
beatSignal1 = beatSignal1 .* exp(-1j * 2 * pi * u * t * tau_system1);
beatSignal2 = beatSignal2 .* exp(-1j * 2 * pi * u * t * tau_system2);

N = 2^24;
f_beat = (-N/2:N/2-1) * (Fs / N);
beatFFT1 = fftshift(fft(beatSignal1, N));
beatFFT1_pos = abs(beatFFT1(N/2+1:end));
beatFFT1_neg = abs(beatFFT1(1:N/2));
[max_val_pos1, idxPos1] = max(beatFFT1_pos);
[max_val_neg1, idxNeg1] = max(beatFFT1_neg);
beatFFT2 = fftshift(fft(beatSignal2, N));
beatFFT2_pos = abs(beatFFT2(N/2+1:end));
beatFFT2_neg = abs(beatFFT2(1:N/2));
[max_val_pos2, idxPos2] = max(beatFFT2_pos);
[max_val_neg2, idxNeg2] = max(beatFFT2_neg);

if max_val_pos1 >= max_val_neg1
    fb1 = abs(f_beat(N/2 + idxPos1));
else
    fb1 = abs(f_beat(idxNeg1));
end

if max_val_pos2 >= max_val_neg2
    fb2 = abs(f_beat(N/2 + idxPos2));
else
    fb2 = abs(f_beat(idxNeg2));
end

%fb = fb2;              % To calibrate check fb1 and fb2 with 1ft cable
fb = (fb1 + fb2) / 2;

tau_from_fb = fb / u;
range_fb = (c * tau_from_fb) / 2;
range_fb_ft = range_fb * 3.28084;
range = (range_fb + range_td)/2;
range_ft = range * 3.28084;

% FFT
Nfft = length(sti(:,1));                   % Number of samples
f_rx = linspace(-Fs/2,Fs/2,Nfft);  % Frequency axis centered -Fs/2 to Fs/2 and shifted
rxFFT1 = fftshift(abs(fft(sri(:,1), Nfft)));
rxFFT2 = fftshift(abs(fft(sri(:,2), Nfft)));

% MUSIC Algorithm for DOA Estimation
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

% Generate Heatmap
thetaScan = -90:0.1:90;
range_axis = linspace(0, 40, 1000);
[ThetaGrid, RangeGrid] = meshgrid(thetaScan, range_axis);
TargetMap = zeros(size(RangeGrid));

kernelSize = 51;
halfSize = floor(kernelSize / 2);
sigma = kernelSize / 9;
[x, y] = meshgrid(-halfSize:halfSize, -halfSize:halfSize);
intensityKernel = exp(-(x.^2 + y.^2) / (2 * sigma^2));
intensityKernel = intensityKernel / max(intensityKernel(:));

% Find location for time-domain estimated range
[~, angleIdx1] = min(abs(thetaScan - AoA)); % Location on x axis
[~, rangeIdx1] = min(abs(range_axis - range_td_ft)); % Location on y axis

% Find location for frequency-domain estimated range
[~, rangeIdx2] = min(abs(range_axis - range_fb_ft));
[~, angleIdx2] = min(abs(thetaScan - AoA)); 

% Average estimated range
[~, rangeIdx3] = min(abs(range_axis - range_ft));
[~, angleIdx3] = min(abs(thetaScan - AoA)); 

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

end
