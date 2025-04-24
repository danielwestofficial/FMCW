function [ret] = save_sc16q11_MIMO(filename, signal)
% SAVE_SC16Q11 Write a normalized complex signal (single or multi-channel)
% to a binary file in the bladeRF "SC16 Q11" format.
%
% Args:
%   filename (str): Target filename. File will be overwritten if it exists.
%   signal (complex matrix): Complex signal with each column representing
%       a channel. Values must be in [-1.0, 1.0).
%
% Returns:
%   ret (int): 1 on success, 0 on failure.

    ret = 0;

    [num_samples, num_channels] = size(signal);
    if num_channels < 1
        error('Input signal must have at least one channel.');
    end

    [f, err_msg] = fopen(filename, 'w', 'ieee-le');
    if f == -1
        error('Failed to open file "%s": %s', filename, err_msg);
    end

    try
        
        sig_i = round(real(signal) .* 2048.0);
        sig_i(sig_i > 2047)  = 2047; % Clamped to prevent overflow
        sig_i(sig_i < -2048) = -2048; % Clamped to prevent overflow
        
        sig_q = round(imag(signal) .* 2048.0);
        sig_q(sig_q > 2047)  = 2047; % Real and imaginary parts are interleaved
        sig_q(sig_q < -2048) = -2048; % Real and imaginary parts are interleaved

        % Interleave I and Q components for all channels
        sig_out = zeros(2 * num_samples * num_channels, 1, 'int16');

        % Interleave indexing 
        for n = 1:num_samples
            for ch = 1:num_channels
                base_idx = (n - 1) * (2 * num_channels) + (ch - 1) * 2 + 1;
                sig_out(base_idx)     = sig_i(n, ch); % I component
                sig_out(base_idx + 1) = sig_q(n, ch); % Q component
            end
        end

        count = fwrite(f, sig_out, 'int16');
        fclose(f);

        if count == length(sig_out)
            ret = 1; % Success
        else
            warning('Failed to write all data to file "%s".', filename);
        end
    catch ME
        fclose(f); 
        rethrow(ME);
    end
end