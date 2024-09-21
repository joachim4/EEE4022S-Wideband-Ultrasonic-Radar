% Parameters
Fs = 192000; % Sampling frequency (192 kHz)
f = 40000;   % Frequency of sine wave (40 kHz)
duration = 2; % Duration in seconds

% Time vector
t = 0:1/Fs:duration;

% Generate 40 kHz sine wave
y = 0.5 * sin(2 * pi * f * t);

% Play the sound
sound(y, Fs);
