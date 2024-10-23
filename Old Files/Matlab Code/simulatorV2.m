%% Clear variables and command window
clear all;
close all;
clc;

%% Define constants and parameters
SpeedSoundWave_ms = 343;             % [m/s]  -> Speed of sound wave
Fc_Hz = 40000;                      % [Hz]   -> Carrier frequency
TimeDuration_s = 10;                % [s]    -> Duration of signal transmission and reception
Fs = 120000;                        % [Hz]   -> Sampling rate

%% Generate the transmit signal
Ts = 1/Fs;                          % Sampling period
t = 0:Ts:(TimeDuration_s);          % Time vector for pulse

% Pure sinusoid for transmission
TxSignal = sin(2 * pi * Fc_Hz * t);

%% Play out transmit signal through the speakers
soundsc(TxSignal, Fs, 24);  % Transmit the signal

%% Record received samples from the microphone
RecLength_samples = length(TxSignal);
RecLength_s = RecLength_samples * Ts; 
recObj = audiorecorder(Fs, 24, 1);
recordblocking(recObj, RecLength_s);  % Record audio for a fixed number of seconds
RX_signal = getaudiodata(recObj);     % Store recorded audio signal in double-precision array

%% Plot transmit and received signals
figure; axes('fontsize', 12);
subplot(2,1,1);
plot(t, TxSignal); % Plot transmit signal
xlabel('Time (s)', 'fontsize', 12);
ylabel('Amplitude (linear)', 'fontsize', 12);
title('Transmit Signal', 'fontsize', 12);
grid on;

subplot(2,1,2);
plot(t, RX_signal); % Plot received signal
xlabel('Time (s)', 'fontsize', 12);
ylabel('Amplitude (linear)', 'fontsize', 12);
title('Received Signal', 'fontsize', 12);
grid on;

% Downmix to baseband for I and Q 
cos_local_osc = cos(2 * pi * Fc_Hz * t); 
sin_local_osc = sin(2 * pi * Fc_Hz * t);  

% Multiply received signal with cosine and sine to get I and Q components
I_signal = RX_signal .* cos_local_osc';  % I component 
Q_signal = RX_signal .* sin_local_osc';  % Q component 

%% Low-pass filter design
fcutoff = 20000;  % Cutoff frequency at 20 kHz
[b, a] = butter(5, fcutoff / (Fs / 2), 'low');  % 5th order Butterworth filter

% Low-pass filter 
I_baseband = filter(b, a, I_signal);
Q_baseband = filter(b, a, Q_signal);

% Combine components to get signal at baseband
baseband_signal = I_baseband + (1i * Q_baseband);

%% Spectrogram to analyze Doppler frequencies
window_length = 1024;  % Length of each segment
overlap = 512;         % Overlap between segments
nfft = 2048;           % Number of FFT points

figure;
spectrogram(baseband_signal, window_length, overlap, nfft, Fs, 'yaxis');
title('Spectrogram of Received Signal (Baseband I/Q)');
ylabel('Frequency (Hz)');
xlabel('Time (s)');
colorbar;


% Parameters for spectrogram
W = 1024; % Frame length - 1024 samples
O = 0.5 * W; % 50% overlap
nfft = 1024; % FFT length
w = hamming(1, W);                                                     

[S, TimeAxis_s, FrequencyAxis_Hz] = JoachimSpectrogram(y, w, W, O, fs);

% Spectrogram plot
imagesc(TimeAxis_s, FrequencyAxis_Hz, 20*log10(abs(S)));              
colorbar;
xlabel('Time');
ylabel('Frequency');
axis xy;

function [S, TimeAxis_s, FrequencyAxis_Hz] = JoachimSpectrogram(x, w, W, O, fs)

% Number of frames
N_f = ((length(x)-W)/O)+1;

FrequencyAxis_Hz = (-W/2:1:(W/2-1))*fs/W;      
Time_s = (0:1:length(x))*1/fs;              

S = zeros(W, N_f);
TimeAxis_s = zeros(1, N_f);                    

% Computing spectrogram
    for k = 1:N_f
        StartIdx = (k - 1) * (W - O) + 1;
        StopIdx = StartIdx + W - 1;
        frame = x(StartIdx : StopIdx);
        frame = frame .* w;
        fft_frame = fftshift(fft(frame, W));
        % S(:, k) = abs(fft_frame).^2;
        S(:, k) = fft_frame;                % Dr Abdul Gaffar. Only use abs( ) when plotting. Keep the complex values here.
        TimeAxis_s(k) = (Time_s(StartIdx) + Time_s(StopIdx)) / 2;
    end
end
