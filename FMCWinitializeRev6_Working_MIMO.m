clear all
close all
clc

%% Step 1: Chirp Parameters
% Values used to define the generated chirp waveform

Fs = 1e6; % Sample Frequency works down to 1.2 MHz but most accurate at 2 MHz (lower sampling saves memory)
frequencyStart = 2.425e9; % Frequency at beginning of chirp
frequencyEnd = 2.475e9; % Frequency at end of chirp
BW = (frequencyEnd - frequencyStart); % Bandwidth of 50 Mhz allows extraction of beat frequency
%BW = Fs/2; % Lowering BW fulfills Nyquist Shannon Theorem
fo = (frequencyStart + frequencyEnd)/2; % Center Frequency
tau = 1e-3; % Sweep time of a single chirp
u = (frequencyEnd-frequencyStart)/tau; % instantaneous slope
A = 0.5; % Amplitude of signal
phi = 0; % Phase offset (used for innacuracies in phase between channels)
w = 2*pi*frequencyStart; % omega (normalized frequency)

%% Step 2: Generate Waveform
% This section provides the time vectors, instantaneous phase, and complex
% chirp signal. 

% Sample Sizes
%NUM_SECONDS = 0.1;               % Duration in seconds
num_samples = 1024;               % Number of samples to receive per buffer

% Time vectors
tTX = (0:1/num_samples:1-1/num_samples); % Time vector
%tTX = 0 : (1/Fs) : (NUM_SECONDS - 1/Fs); % Time vector for all samples

% Time vector for single chirp
tTXchirp = 0:1/Fs:tau-1/Fs;  
tTXchirp = repmat(tTXchirp, 1, num_samples);
tTXchirp = tTXchirp(1:length(tTX)); % Truncate to match tTX

% Instantaneous frequency of the chirp 
thetaTX = 2*pi*(frequencyStart * tTXchirp + 0.5 * u * tTXchirp.^2); 
thetaTX = repmat(thetaTX, 1, num_samples);
thetaTX = thetaTX(1:length(tTX)); % Truncate to match tTX

% Generate complex chirp
st = A * exp(1i * (w*tTXchirp + thetaTX + phi)); % Complex chirp


%% Step 3: Initialize bladeRF Device for 2 TX and 2 RX Channels (MIMO)

try 
    % Create a bladeRF device handle. Use the serial number if needed (e.g., '*:serial=43b').
    device = bladeRF_MIMO('*:serial=608');
    
    % Print device information
    disp('Device initialized successfully.');
    fprintf('Device Serial: %s\n', device.info.serial);
    fprintf('FPGA Size: %s\n', device.info.fpga_size);
    fprintf('USB Speed: %s\n', device.info.usb_speed);
    
    % Configure TX and RX Channels
    device.tx.channel = 'TX'; % Enable TX1 and TX2
    device.rx.channel = 'RX'; % Enable RX1 and RX2
    fprintf('Active TX Channels: %s\n', strjoin(device.tx.channel, ', '));
    fprintf('Active RX Channels: %s\n', strjoin(device.rx.channel, ', '));

    % Initialize PA and LNA
    %device.misc.attachExpansion('xb100'); % BT100 for TX PA
    %device.misc.attachExpansion('xb200'); % BT200 for RX LNA

    device.loopback = 'NONE'; % Set loopback mode (NONE, FIRMWARE, RFIC_BIST)

    %%%%%%%%%%%% The RX gain requires all channels to be formatted like this %%%%%%%%%%%%  
    % RX Configuration
    device.rx.frequency = fo;  % Center frequency
    device.rx.samplerate = Fs;   % Sample rate
    device.rx.bandwidth = BW;     % Bandwidth
    device.rx.agc = 'manual';
    device.rx.gain = [15, 15];     % Gain for both channels

    % TX Configuration
    device.tx.frequency = fo;  % Center frequency
    device.tx.samplerate = Fs;   % Sample rate
    device.tx.bandwidth = BW;     % Bandwidth
    device.tx.gain = [20, 20];     % Gain for both channels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('Default TX Gain: [%d, %d]\n', device.tx.gain(1), device.tx.gain(2));
    fprintf('Default RX Gain: [%d, %d]\n', device.rx.gain(1), device.rx.gain(2));

    % Print RX and TX configuration
    fprintf('RX Config: Freq=%.2f GHz, Rate=%.2f MSps, BW=%.2f MHz, Gain=[%d,%d]\n', ...
        device.rx.frequency / 1e9, device.rx.samplerate / 1e6, device.rx.bandwidth / 1e6, ...
        device.rx.gain(1), device.rx.gain(2));
    fprintf('TX Config: Freq=%.2f GHz, Rate=%.2f MSps, BW=%.2f MHz, Gain=[%d,%d]\n', ...
        device.tx.frequency / 1e9, device.tx.samplerate / 1e6, device.tx.bandwidth / 1e6, ...
        device.tx.gain(1), device.tx.gain(2));
catch ME
    disp('Error during device initialization:');
    disp(getReport(ME, 'extended'));
    %disp(ME.message);
    return;
end


%% Step 4: Configure RX and TX Streams (WORKING)

num_buffers = 16;               % Number of buffers
buffer_size = num_samples;      % Buffer size in samples
num_transfers = 8;              % Number of active USB transfers
timeout = 5000;                 % Timeout in milliseconds

% Configure RX streaming
device.rx.config.num_buffers = num_buffers;
device.rx.config.buffer_size = buffer_size;
device.rx.config.num_transfers = num_transfers;
device.rx.config.timeout_ms = timeout;

% Configure TX streaming
device.tx.config.num_buffers = num_buffers;
device.tx.config.buffer_size = buffer_size;
device.tx.config.num_transfers = num_transfers;
device.tx.config.timeout_ms = timeout;

% Configure the RX buffer
sri = zeros(num_samples, 2); 
rx_channel0 = sri(:, 1); 
rx_channel1 = sri(:, 2); 
sri = [rx_channel0(:), rx_channel1(:)]; 

%% Start RX and TX Streams

% Configure the TX buffer
sti = zeros(num_samples, 2); % Preallocate TX buffer
tx_channel0 = A * exp(1i * (w*tTXchirp + thetaTX + phi))';
%tx_channel1 = A * exp(1i * (w*tTXchirp + thetaTX + phi + pi)); % tx1 will reflect phase shift from antenna array
tx_channel1 = A * exp(1i * (w*tTXchirp + thetaTX + phi))';

% Combine into [N, 2] matrix
sti = [tx_channel0(:), tx_channel1(:)]; % Ensure two columns
max_amp = max(abs(sti(:)));
sti = sti / max_amp; % Normalize to avoid clipping

% Start Streams 
device.rx.start();
pause(0.01); 
device.tx.start();

%% Run

try
    for i = 1:64
        tic;
        device.transmit(sti);
        pause(0.0001); % This value is sensitive to acheive samples. Can make it seem like you are not receiving because you missed transmittion!
        sri = device.receive(num_samples);
        elapsed = toc; 
        fprintf('Iteration %d elapsed time: %.4f seconds\n', i, elapsed);
    end
catch ME
    disp('Error during TX/RX operation:');
    %disp(ME.message);
    disp(getReport(ME, 'extended'));
    device.tx.stop();
    device.rx.stop();
    device.reset(); % Reset the device
    return;
end

% Stop TX and RX
device.tx.stop();
device.rx.stop();

% Save Receive Signal as a SC16 Q11 File
%save_sc16q11('C:\Users\danie\OneDrive\Documents\OIT(FALL.2024)\ENGR596\Week10\HardwareInitialize(bladeRF_MIMOrev1)\tmp\sri.sc16q11', sri)


%% Step 5: Verify and Plot Received Data

% Load the samples
%sri = load_sc16q11('sri.sc16q11');

% Check if the received data is valid
if all(sri(:) == 0)
    warning('Received data is all zeros. Check the gain settings and loopback configuration.');
end

% Create a corrected time vector
%tTX = (0 : (length(sti) - 1)) / Fs;
%tRX = (0 : (length(sri) - 1)) / Fs;
tTX = (0 : size(sti, 1) - 1) / Fs;
tRX = (0 : size(sri, 1) - 1) / Fs;

% Plot TX signal for verification
figure;
subplot(2,3,1);
plot(tTX, sti); % Plot and zoom
title('Real Chirp Signal');
xlabel('Time (s)');
xlim([0 1e-4]);
ylabel('Amplitude');

subplot(2,3,2);
plot(tTX, real(sti));
title('Complex Chirp Signal - Real Part');
xlabel('Time (s)');
xlim([0 1e-4]);
ylabel('Amplitude');

subplot(2,3,3);
plot(tTX, imag(sti));
title('Complex Chirp Signal - Imaginary Part');
xlabel('Time (s)');
xlim([0 1e-4]);
ylabel('Amplitude');

% Plot RX signal for verification
subplot(2,3,4);
plot(tRX, sri); % Plot and zoom
title('Real Chirp Signal');
xlabel('Time (s)');
xlim([0 1e-4]);
ylabel('Amplitude');

subplot(2,3,5);
plot(tRX, real(sri));
title('Complex Chirp Signal - Real Part');
xlabel('Time (s)');
xlim([0 1e-4]);
ylabel('Amplitude');

subplot(2,3,6);
plot(tRX, imag(sri));
title('Complex Chirp Signal - Imaginary Part');
xlabel('Time (s)');
xlim([0 1e-4]);
ylabel('Amplitude');

disp('Loopback test complete.');

%% Cleanup
device.tx.stop();
device.rx.stop();
clear device;
disp('Cleanup completed.');
