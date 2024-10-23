%% Clear variables and command window
clear all;
close all;
clc;

%% Load the saved data
load('received_signal.mat');

%% Parameters
chunk_size = 1000; % Process 1000 samples at a time
num_chunks = ceil(length(RX_signal) / chunk_size);

%% Demodulation and filtering in small chunks
fcutoff = 20000;
[b, a] = butter(5, fcutoff / (Fs / 2), 'low');

% Initialize filter states
zi_I = zeros(max(length(a),length(b))-1, 1);
zi_Q = zeros(max(length(a),length(b))-1, 1);

%% Process data in small chunks
window_length = 256;
overlap = 128;
nfft = 256;
freq_range = [0, fcutoff];

% Initialize spectrogram accumulator
S_acc = zeros(nfft/2+1, ceil(length(RX_signal)/window_length));
col_idx = 1;

for i = 1:num_chunks
    start_idx = (i-1)*chunk_size + 1;
    end_idx = min(i*chunk_size, length(RX_signal));
    
    chunk_t = t(start_idx:end_idx);
    chunk_RX = RX_signal(start_idx:end_idx);
    
    cos_component = cos(2 * pi * Fc_Hz * chunk_t);
    sin_component = sin(2 * pi * Fc_Hz * chunk_t);
    
    I_signal = chunk_RX .* cos_component;
    Q_signal = chunk_RX .* sin_component;
    
    [I_baseband, zi_I] = filter(b, a, I_signal, zi_I);
    [Q_baseband, zi_Q] = filter(b, a, Q_signal, zi_Q);
    
    baseband_signal = I_baseband + (1i * Q_baseband);
    
    % Calculate spectrogram for this chunk
    [S, F, T] = spectrogram(baseband_signal, window_length, overlap, nfft, Fs, 'yaxis');
    
    % Accumulate spectrogram
    num_cols = size(S, 2);
    if col_idx + num_cols - 1 > size(S_acc, 2)
        S_acc = [S_acc, zeros(size(S_acc, 1), num_cols)];  % Extend S_acc if needed
    end
    S_acc(:, col_idx:col_idx+num_cols-1) = S_acc(:, col_idx:col_idx+num_cols-1) + abs(S);
    col_idx = col_idx + num_cols;
end

% Trim excess columns from S_acc
S_acc = S_acc(:, 1:col_idx-1);

%% Spectrogram plot
T = (0:size(S_acc,2)-1) * (window_length - overlap) / Fs;
F = linspace(0, Fs/2, nfft/2+1);

figure;
imagesc(T, F, 20*log10(S_acc));
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram of Received Signal');
axis xy;
ylim(freq_range);

% Display processing information
disp(['Processed ', num2str(length(RX_signal)), ' samples in ', num2str(num_chunks), ' chunks.']);
disp(['Resulting spectrogram size: ', num2str(size(S_acc))]);