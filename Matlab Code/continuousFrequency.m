% Parameters
Fs = 44100; % Sampling frequency (Hz)
nBits = 16; % Bit depth
nChannels = 1; % Number of channels (1 for mono, 2 for stereo)
chunkDuration = 0.1; % Duration of each chunk (seconds)

% Create an audiorecorder object
recObj = audiorecorder(Fs, nBits, nChannels);

% Initialize the time-domain plot
figure;
subplot(2,1,1); % Top subplot for time-domain
hPlotTime = plot(NaN, NaN);
title('Live Microphone Input - Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Initialize the frequency-domain plot
subplot(2,1,2); % Bottom subplot for frequency-domain
hPlotFreq = plot(NaN, NaN);
title('Live Microphone Input - Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Start continuous recording
disp('Start speaking... Press Ctrl+C to stop.');
record(recObj);

% Continuous plotting loop
while true
    % Record a chunk of audio
    pause(chunkDuration); % Pause to allow recording of the chunk
    audioData = getaudiodata(recObj);

    % Update the time-domain plot
    t = (0:length(audioData)-1) / Fs; % Time vector
    set(hPlotTime, 'XData', t, 'YData', audioData);
    
    % Calculate and update the frequency-domain plot
    L = length(audioData); 
    Y = fft(audioData); % FFT of the audio data
    P2 = abs(Y/L); % Two-sided spectrum
    P1 = P2(1:L/2+1); % Single-sided spectrum
    f = Fs*(0:(L/2))/L; % Frequency vector
    
    set(hPlotFreq, 'XData', f, 'YData', P1);
    
    % Force MATLAB to update the figure window
    drawnow;
end

% Stop recording (if needed, this line will not be reached)
stop(recObj);
