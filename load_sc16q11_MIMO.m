function [ signal, signal_i, signal_q ] = load_sc16q11_MIMO(filename, num_channels)
% LOAD_SC16Q11_MIMO Load a binary file containing multiple channels of 
% complex samples in the bladeRF "SC16 Q11" format.
%
% Args:
%   filename (str): File to load.
%   num_channels (int): Number of channels in the file.
%
% Returns:
%   signal (complex matrix): Loaded signal where each column is a channel.

    % Open file
    [f, err_msg] = fopen(filename, 'r', 'ieee-le');
    if f == -1
        error("Failed to open file: " + err_msg);
    end

    data = fread(f, 'int16');
    fclose(f);

    % Debugging
    %disp('Raw data (first 16 samples):');
    %disp(data(1:16));

    total_samples = length(data);
    if mod(total_samples, 2 * num_channels) ~= 0
        error('File length is not compatible with the expected number of channels.');
    end

    num_samples = total_samples / (2 * num_channels); % Compute number of IQ samples per channel

    signal = zeros(num_samples, num_channels); % Initialize output matrix

    % Extract I/Q data per channel
    for ch = 1:num_channels
        i_idx = (ch - 1) * 2 + 1 : 2 * num_channels : total_samples;
        q_idx = i_idx + 1;
        
        signal_i = double(data(i_idx)) / 2048;
        signal_q = double(data(q_idx)) / 2048;
        
        signal(:, ch) = signal_i + 1j * signal_q;
    end
end