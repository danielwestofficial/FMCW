%% FMCW Script
clear all;
close all;
clc;

% Radar Parameters
Fs = 60e6;
Ts = 1/Fs;
BW = 12e6;
frequencyStart = 2.45e9 - BW/2;
frequencyEnd = 2.45e9 + BW/2;
fo = (frequencyStart + frequencyEnd)/2;
tau =10e-6;
u = BW / tau;
A = 1;
phi = 0;

NUM_REPEATS = 3;
NUM_SECONDS = tau * NUM_REPEATS;
NUM_SAMPLES = NUM_SECONDS * Fs;

% Initializing Chirp
tchirp = 0:1/Fs:tau - 1/Fs;
thetaTX = 2*pi*(tchirp + 0.5 * u * tchirp.^2);
tx_channel0 = A * exp(1i * (thetaTX + phi));
tx_channel1 = tx_channel0;
tx_channel0_repeated = repmat(tx_channel0, 1, NUM_REPEATS); 
tx_channel1_repeated = repmat(tx_channel1, 1, NUM_REPEATS); 
sti = [tx_channel0_repeated(:), tx_channel1_repeated(:)];   % Preallocate buffer

% Apply function to save signal to sc16q11 file
save_sc16q11_MIMO('R:\Temp\chirpLoop.sc16q11', sti); 

% Update time vector for plotting
t = linspace(0, NUM_SECONDS, length(tx_channel0_repeated));

% Plot Real and Imaginary Parts of Chirp Signal
figure;
subplot(2,1,1);
plot(t, real(sti(:,1)), 'b', t, imag(sti(:,1)), 'r');
title('Time-Domain LFM Upchirp Signal');
xlabel('Time (s)');
xlim([0 30e-6]);
ylabel('Amplitude');
legend('Real Part', 'Imaginary Part');
grid on;

subplot(2,1,2);
plot(t, real(sti(:,2)), 'b', t, imag(sti(:,2)), 'r');
title('Time-Domain LFM Upchirp Signal');
xlabel('Time (s)');
xlim([0 30e-6]);
ylabel('Amplitude');
legend('Real Part', 'Imaginary Part');
grid on;

%% Run Master & Slave Windows Batch Files
% Further Directions in bladeRF CLI
system('runMaster.bat');
system('runSlave.bat');
