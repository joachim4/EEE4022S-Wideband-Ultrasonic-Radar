% Parameters
Fs = 192000; % Sampling frequency (Hz)
nBits = 24; % Bit depth
nChannels = 2; % Number of channels (1 for mono, 2 for stereo)
duration = 5; % Recording duration (seconds)

% Create an audiorecorder object
recObj = audiorecorder(Fs, nBits, nChannels);

% Start recording
disp('Start speaking...');
recordblocking(recObj, duration);
disp('Recording complete.');

% Retrieve the recorded audio data
audioData = getaudiodata(recObj);

% Plot the recorded audio signal
t = (0:length(audioData)-1) / Fs; % Time vector
figure;
plot(t, audioData);
title('Recorded Audio Signal');
xlabel('Time (s)');
ylabel('Amplitude');
