% Parameters
Fs = 192000; % Sampling frequency (Hz)
nBits = 24; % Bit depth
nChannels = 1; % Number of channels (1 for mono, 2 for stereo)
chunkDuration = 0.1; % Duration of each chunk (seconds)

% Create an audiorecorder object
recObj = audiorecorder(Fs, nBits, nChannels);

% Initialize the plot
figure;
hPlot = plot(NaN, NaN);
title('Live Microphone Input');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Start continuous recording
disp('Start speaking... Press Ctrl+C to stop.');
record(recObj);

% Continuous plotting loop
while true
    % Record a chunk of audio
    pause(chunkDuration); % Pause to allow recording of the chunk
    audioData = getaudiodata(recObj);

    % Update the plot
    t = (0:length(audioData)-1) / Fs; % Time vector
    set(hPlot, 'XData', t, 'YData', audioData);
    drawnow;
end

% Stop recording (if needed, this line will not be reached)
stop(recObj);
