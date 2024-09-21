% Parameters
Fs = 192000;         % Sampling frequency (192 kHz is common for high-quality audio cards)
t = 0:1/Fs:1;        % Time vector for 1 second duration
f = 40000;           % Frequency of sine wave (40 kHz)

% Generate sine wave
y = 0.5 * sin(2 * pi * f * t);  % Amplitude of 0.5 to avoid clipping

% Play the sound
player = audioplayer(y, Fs);
play(player);
