%% Clear variables and command window
clear all;
close all;
clc;

%% Define constants and parameters
SpeedSoundWave_ms = 343; % [m/s] -> Speed of sound wave
Fc_Hz = 40000; % [Hz] -> Carrier frequency
TimeDuration_s = 10; % [s] -> Duration of signal transmission and reception
Fs = 120000; % [Hz] -> Sampling rate
v_target_kmh = 60; % [km/h] -> Velocity of the target (for simulation purposes)
v_target_ms = v_target_kmh / 3.6; % [m/s] -> Convert velocity to m/s

%% Calculate Doppler Frequency (for simulation purposes)
f_doppler = (2 * v_target_ms * Fc_Hz) / SpeedSoundWave_ms;
disp(['Actual Doppler frequency: ' num2str(f_doppler) ' Hz']);

%% Generate the transmit signal
Ts = 1/Fs; % Sampling period
t = 0:Ts:(TimeDuration_s);
TxSignal = sin(2 * pi * Fc_Hz * t);

%% Simulate the received signal with Doppler effect
%RX_signal = sin(2 * pi * (Fc_Hz + f_doppler) * t); % Doppler-shifted received signal

%% Play out transmit signal through the speakers
soundsc(TxSignal, Fs, 24);  % Transmit the signal

%% Record received samples from the microphone
RecLength_samples = length(TxSignal);
RecLength_s = RecLength_samples * Ts; 
recObj = audiorecorder(Fs, 24, 1);
recordblocking(recObj, RecLength_s);  % Record audio for a fixed number of seconds
RX_signal = getaudiodata(recObj).';     % Store recorded audio signal in double-precision array
 
clear recObj;

%% Plot transmit and received signals
figure; axes('fontsize', 12);
subplot(2,1,1);
plot(t, TxSignal);
xlabel('Time (s)', 'fontsize', 12);
ylabel('Amplitude (linear)', 'fontsize', 12);
title('Transmit Signal', 'fontsize', 12);
grid on;
subplot(2,1,2);
plot(t, RX_signal);
xlabel('Time (s)', 'fontsize', 12);
ylabel('Amplitude (linear)', 'fontsize', 12);
title('Received Signal (with Doppler Shift)', 'fontsize', 12);
grid on;

%% Quadrature Demodulation
cos_component = cos(2 * pi * Fc_Hz * t);
sin_component = sin(2 * pi * Fc_Hz * t);
I_signal = RX_signal .* cos_component;
Q_signal = RX_signal .* sin_component;

%% Low-pass filter design
fcutoff = 20000;
[b, a] = butter(5, fcutoff / (Fs / 2), 'low');

% Low-pass filter
I_baseband = filter(b, a, I_signal);
Q_baseband = filter(b, a, Q_signal);
baseband_signal = I_baseband + (1i * Q_baseband);

%% Spectrogram analysis
% Parameters for spectrogram
W = 1024; % Frame length - 1024 samples
O = 0.5 * W; % 50% overlap
nfft = 1024; % FFT length
w = hamming(W);

[S, TimeAxis_s, FrequencyAxis_Hz] = JoachimSpectrogram(baseband_signal, w, W, O, Fs, nfft);

% Spectrogram plot
figure;
imagesc(TimeAxis_s, FrequencyAxis_Hz, 20*log10(abs(S)));
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram of Received Signal');
axis xy;

%% Speed Detection Algorithm
% 1. Find the peak frequency in each time frame
[~, peak_indices] = max(abs(S));
peak_frequencies = FrequencyAxis_Hz(peak_indices);

% 2. Calculate the median peak frequency to reduce noise effects
median_peak_frequency = median(peak_frequencies);

% 3. Calculate the estimated speed using the Doppler equation
estimated_speed_ms = (median_peak_frequency * SpeedSoundWave_ms) / (2 * Fc_Hz);
estimated_speed_kmh = estimated_speed_ms * 3.6;

disp(['Estimated target speed: ' num2str(estimated_speed_kmh) ' km/h']);
disp(['Actual target speed: ' num2str(v_target_kmh) ' km/h']);
disp(['Error: ' num2str(abs(abs(estimated_speed_kmh) - abs(v_target_kmh))) ' km/h']);

%% Spectrogram Function
function [S, TimeAxis_s, FrequencyAxis_Hz] = JoachimSpectrogram(x, w, W, O, fs, nfft)
% Ensure x is a column vector
 x = x(:);
% Number of frames
 N_f = floor((length(x)-W)/O)+1;
% Pre-allocate the spectrogram matrix
 S = zeros(nfft, N_f);
% Time axis
 TimeAxis_s = ((0:N_f-1) * O + W/2) / fs;
% Frequency axis
if mod(nfft, 2) == 0
 FrequencyAxis_Hz = (-nfft/2:nfft/2-1) * fs / nfft;
else
 FrequencyAxis_Hz = (-(nfft-1)/2:(nfft-1)/2) * fs / nfft;
end
% Computing spectrogram
for k = 1:N_f
 StartIdx = (k - 1) * O + 1;
 StopIdx = StartIdx + W - 1;
 frame = x(StartIdx : StopIdx);
 frame = frame .* w;
 fft_frame = fftshift(fft(frame, nfft));
 S(:, k) = fft_frame;
end
end