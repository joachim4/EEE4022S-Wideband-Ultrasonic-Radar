2% Parameters
fs = 120e3;            % Sampling frequency (120 kHz)
T = 1;                 % Duration of the signal in seconds
fc = 40e3;             % Carrier frequency (40 kHz)
t = (0:1/fs:T-1/fs)';  % Time vector

% Transmitted signal
tx_signal = cos(2 * pi * fc * t);

% Doppler Effect Simulation 
v_objects = [5, 10, 15];   % Velocities of objects (m/s)
fd = 2 * fc * v_objects / 343;  % Doppler shifts (Hz), assuming speed of sound 343 m/s

% Create the received signal with Doppler shifts
rx_signal = zeros(size(tx_signal));
for i = 1:length(v_objects)
    rx_signal = rx_signal + cos(2 * pi * (fc + fd(i)) * t);  % Sum of received signals
end

% Downmix to baseband for I (in-phase) and Q (quadrature-phase) components
cos_local_osc = cos(2 * pi * fc * t);  % Local oscillator - cosine
sin_local_osc = sin(2 * pi * fc * t);  % Local oscillator - sine

I_signal = rx_signal .* cos_local_osc;  % Multiply received signal with cosine (I component)
Q_signal = rx_signal .* sin_local_osc;  % Multiply received signal with sine (Q component)

% Low-pass filter design
fcutoff = 20e3;  % Cutoff frequency (20 kHz, covering the Doppler range)
[b, a] = butter(5, fcutoff / (fs / 2), 'low');  % 5th order Butterworth filter

% Apply low-pass filter to extract baseband signals
I_baseband = filter(b, a, I_signal);
Q_baseband = filter(b, a, Q_signal);

% Combine I and Q components to form the complex baseband signal
baseband_signal = I_baseband + 1i * Q_baseband;

% Spectrogram to analyze Doppler frequencies
window_length = 1024;  % Length of each segment
overlap = 512;         % Overlap between segments
nfft = 2048;           % Number of FFT points

figure;
spectrogram(baseband_signal, window_length, overlap, nfft, fs, 'yaxis');
title('Spectrogram of Received Signal (Baseband I/Q)');
ylabel('Frequency (Hz)');
xlabel('Time (s)');
colorbar;
