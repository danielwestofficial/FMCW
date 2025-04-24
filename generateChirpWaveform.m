function [sti,NUM_REPEATS] = generateChirpWaveform(transmitFile, Fs, BW, tau, NUM_REPEATS)
% generateChirpWaveform generates a baseband FMCW chirp waveform and saves
% it to an SC16Q11 file for dual-channel transmission.
%
% INPUTS:
%   transmitFile    - Full file path to save the SC16Q11 waveform (string)
%   Fs            - Sampling rate (Hz), default = 60e6
%   BW            - Bandwidth (Hz), default = 12e6
%   tau           - Chirp duration (sec), default = 10e-6
%   NUM_REPEATS   - Number of chirp repetitions, default = 3
%
% OUTPUT:
%   sti           - [N x 2] complex baseband chirp matrix (MIMO format)

    if nargin < 5, NUM_REPEATS = 1; end
    if nargin < 4, tau = 10e-6; end
    if nargin < 3, BW = 30e6; end
    if nargin < 2, Fs = 60e6; end
    if nargin < 1, transmitFile = 'R:\Temp\transmit.sc16q11'; end

    % Derived Parameters
    Ts = 1/Fs;
    frequencyStart = 2.45e9 - BW/2;
    frequencyEnd = 2.45e9 + BW/2;
    fo = (frequencyStart + frequencyEnd) / 2;
    u = BW / tau;     % Chirp slope
    A = 1;
    phi = 0;

    % Time Vectors and Chirp Signal
    NUM_SECONDS = tau * NUM_REPEATS;
    tchirp = 0:Ts:tau - Ts;
    thetaTX = 2*pi*(tchirp + 0.5 * u * tchirp.^2);
    tx_chirp = A * exp(1i * (thetaTX + phi));

    % Repeat for NUM_REPEATS and format for MIMO
    tx_channel0_repeated = repmat(tx_chirp, 1, NUM_REPEATS);
    tx_channel1_repeated = tx_channel0_repeated;
    sti = [tx_channel0_repeated(:), tx_channel1_repeated(:)];

    % Save to SC16Q11 File
    save_sc16q11_MIMO(transmitFile, sti);

    % Optional Plot (for debugging or verification)
    t = linspace(0, NUM_SECONDS, length(tx_channel0_repeated));
    figure;
    subplot(2,1,1);
    plot(t, real(sti(:,1)), 'b', t, imag(sti(:,1)), 'r');
    title('TX Channel 0 - Time-Domain LFM Upchirp');
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([0 10e-6]);
    legend('Real', 'Imag');
    grid on;

    subplot(2,1,2);
    plot(t, real(sti(:,2)), 'b', t, imag(sti(:,2)), 'r');
    title('TX Channel 1 - Time-Domain LFM Upchirp');
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([0 10e-6]);
    legend('Real', 'Imag');
    grid on;
end