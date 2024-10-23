%% This code was created by Joachim Gengan GNGJOA003 - in partial fulfulment of EEE4022S research project
%% The code has been adapted from example code provided by Dr M. Y. Abdul Gaffar

%% Contains code for both simulating a received signal or using the microphone of the computer to receive real signal from sound card
%% Transmit signal is played through the speakers of the computer and filtering is done to the received signal before spectrograms are plotted

clear all;
close all;
clc;

%% Define constants and parameters of the ultrasonic system
SpeedSoundWave_ms = 343;             % [m/s]  -> Speed of sound wave
Fc_Hz = 40000;                       % [Hz]   -> Carrier frequency
TimeDuration_s = 10;                  % [s]    -> Duration of signal transmission and reception
Fs = 192000;                         % [Hz]   -> Sampling rate
v_target_kmh = 36;                   % [km/h] -> Velocity of the target 
v_target_ms = v_target_kmh / 3.6;    % [m/s]  -> Convert velocity to m/s
lamda = SpeedSoundWave_ms/Fc_Hz;

% Creates folders for storing the recorded data
if ~exist('Figures', 'dir')
    mkdir('Figures');
end
if ~exist('Figures\TransmitReceiveSignals', 'dir')
    mkdir('Figures\TransmitReceiveSignals');
end
if ~exist('Figures\Spectrograms', 'dir')
    mkdir('Figures\Spectrograms');
end
if ~exist('Workspaces', 'dir')
    mkdir('Workspaces');
end

%% Unique identifier for this run using current timestamp for easily identifying recorded datasets
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
run_id = getNextRunID();

%% Calculate Doppler Frequency - used for Radar Simulator 
%f_doppler = (2 * v_target_ms * Fc_Hz) / SpeedSoundWave_ms; 
%disp(['Calculated Doppler frequency: ' num2str(f_doppler) ' Hz']);

%% Generate the transmit signal
Ts = 1/Fs;                           
t = 0:Ts:(TimeDuration_s);         
TxSignal = sin(2 * pi * Fc_Hz * t);
soundsc(TxSignal, Fs, 24);  % Transmit the signal through the speakers on computer - connected to sound card

%% Simulate the received signal with Doppler effect (much weaker observed echo) - used for Radar Simulator
%RX_signal = sin(2 * pi * Fc_Hz * t) + 0.01*sin(2 * pi * (Fc_Hz + f_doppler) * t);

%% Record received samples from the microphone
RecLength_samples = length(TxSignal);
RecLength_s = RecLength_samples * Ts; 
recObj = audiorecorder(Fs, 24, 1);
recordblocking(recObj, RecLength_s);  
RX_signal = getaudiodata(recObj).';     % Store recorded audio signal in an array to be used as received signal during further processing
clear recObj; % Make code more efficient, not using object anymore so can free memory

%% Digital bandpass filter
lowFreq = 20000;  % Lower cutoff frequency same as expected for analogue circuitry
highFreq = 60000; % Upper cutoff frequency
[b, a] = butter(6, [lowFreq highFreq]/(Fs/2), 'bandpass');
RX_signal_filtered = filtfilt(b, a, RX_signal);

%% Plot of the transmit and received signals in time domain
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

%% Save the figure as a picture and figure file
saveas(gcf, ['Figures9\TransmitReceiveSignals\TransmitReceiveSignals_' timestamp '_' num2str(run_id) '.png']);
savefig(['Figures9\TransmitReceiveSignals_' timestamp '_' num2str(run_id) '.fig']);

%% FFT and FFTShift
% Compute the FFT of the signals
N = length(t);
f = (-N/2:N/2-1)*(Fs/N);  % Frequency vector 

TxSignal_FFT = fft(TxSignal);
TX_signal_FFTShifted = fftshift(TxSignal_FFT);

RX_signal_FFT = fft(RX_signal_filtered);
RX_signal_FFTShifted = fftshift(RX_signal_FFT);

% FFT shifted signals of Rx and Tx
figure;
subplot(2,1,1);
plot(f, abs(TX_signal_FFTShifted)/N);
xlabel('Frequency (Hz)', 'fontsize', 12);
ylabel('Magnitude', 'fontsize', 12);
title('FFT Shift of Transmit Signal', 'fontsize', 12);
grid on;

subplot(2,1,2);
plot(f, abs(RX_signal_FFTShifted)/N);
xlabel('Frequency (Hz)', 'fontsize', 12);
ylabel('Magnitude', 'fontsize', 12);
title('FFT Shift of Received Signal', 'fontsize', 12);
grid on;

%% Downmixing components needed
cos_component = cos(2 * pi * Fc_Hz * t); 
sin_component = -sin(2 * pi * Fc_Hz * t);  

% Separate signal into In-phase and Quadrature components
I_signal = RX_signal .* cos_component;
Q_signal = RX_signal .* sin_component;  

%% Low-pass filter design
fcutoff = 20000;  
[b, a] = butter(5, fcutoff / (Fs / 2), 'low'); 
[h, f] = freqz(b, a, 1024, Fs);

% Low pass filter frequency response
figure;
plot(f, 20*log10(abs(h)));
grid on;
title('Frequency Response of the Butterworth Filter');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
xlim([0 Fs/2]);

% Low pass filter applied at baseband to attenuate high frequency copies
I_baseband = filter(b, a, I_signal);
Q_baseband = filter(b, a, Q_signal);

baseband_signal = I_baseband + (1i * Q_baseband); 

N_baseband = length(baseband_signal);  
f_baseband = (-N_baseband/2:N_baseband/2-1)*(Fs/N_baseband);  

baseband_FFT = fft(baseband_signal);
baseband_FFTShifted = fftshift(baseband_FFT);

%figure;
%plot(f_baseband, abs(baseband_FFTShifted)/N_baseband);
%xlabel('Frequency (Hz)', 'fontsize', 12);
%ylabel('Magnitude', 'fontsize', 12);
%title('FFT Shift of Baseband Signal', 'fontsize', 12);
%grid on;

% Parameters for spectrogram
W = 1024; % Frame length 
O = 0.5 * W; % 50% overlap
nfft = 1024; % FFT length
w = hamming(W);                                                     

[S, TimeAxis_s, FrequencyAxis_Hz] = JoachimSpectrogram(baseband_signal, w, W, O, Fs, nfft);

% Spectrogram plot for Doppler frequency
figure;
clims = [-50 0];
imagesc(TimeAxis_s, FrequencyAxis_Hz/1e3, 20*log10(abs(S)),clims);              
colorbar;
colormap('jet');
xlabel('Time (s)');
ylabel('Frequency (kHz)');
title('Spectrogram of Received Signal');
axis xy;

saveas(gcf, ['Figures\Spectrograms\Spectrogram_' timestamp '_' num2str(run_id) '.png']);
savefig(['Figures\Spectrograms\Spectrogram_' timestamp '_' num2str(run_id) '.fig']);

%% Plot Zoomed-in spectrogram and convert frequency to velocity (km/hr)
maxSpeed_km_hr = 80;
minSpeed_km_per_hr = -80; 

speed_m_per_sec = FrequencyAxis_Hz*lamda/2;
speed_km_per_hr = speed_m_per_sec*(60*60/1000);
speed_km_per_hr_Idx = find((speed_km_per_hr <= maxSpeed_km_hr) & (speed_km_per_hr >= minSpeed_km_per_hr));

SpeedVectorOfInterest = speed_km_per_hr(speed_km_per_hr_Idx);
S_OfInterest = S(speed_km_per_hr_Idx, :);

S_OfInterestToPlot = abs(S_OfInterest)/max(max(abs(S_OfInterest))); % normalise plot

% Plots of the spectrogram 
clims = [-50 0];
figure; imagesc(t,SpeedVectorOfInterest,20*log10(S_OfInterestToPlot), clims);
xlabel('Time (s)');
ylabel('Speed (km/hr)');
title('Spectrogram: zoomed in (km/hr)');
grid on;
colorbar;
colormap('jet');
axis xy;

saveas(gcf, ['Figures\Spectrograms\SpectrogramKMPH_' timestamp '_' num2str(run_id) '.png']);
savefig(['Figures\Spectrograms\SpectrogramKMPH_' timestamp '_' num2str(run_id) '.fig']);

%% Plot Zoomed-in spectrogram and convert frequency to velocity (m/s)
maxSpeed_m_s = 20;
minSpeed_m_s = -20; 

speed_m_per_sec = FrequencyAxis_Hz*lamda/2;
speed_m_per_sec_Idx = find((speed_m_per_sec <= maxSpeed_m_s) & (speed_m_per_sec >= minSpeed_m_s));

SpeedVectorOfInterest_m_s = speed_m_per_sec(speed_km_per_hr_Idx);
S_OfInterest_m_s = S(speed_m_per_sec_Idx, :);

S_OfInterestToPlot_m_s = abs(S_OfInterest_m_s)/max(max(abs(S_OfInterest_m_s))); % normalise plot

% Plot the spectrogram 
clims = [-50 0];
figure; imagesc(t,SpeedVectorOfInterest_m_s,20*log10(S_OfInterestToPlot_m_s), clims);
xlabel('Time (s)');
ylabel('Speed (m/s)');
title('Spectrogram: zoomed in (m/s)');
grid on;
colorbar;
colormap('jet');
axis xy;

saveas(gcf, ['Figures\Spectrograms\SpectrogramMS_' timestamp '_' num2str(run_id) '.png']);
savefig(['Figures\Spectrograms\SpectrogramKMS_' timestamp '_' num2str(run_id) '.fig']);

save(['Workspaces\SpectrogramWorkspace_' timestamp '_' num2str(run_id) '.mat']);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          

function [S, TimeAxis_s, FrequencyAxis_Hz] = JoachimSpectrogram(x, w, W, O, fs, nfft)
    x = x(:);
    
    N_f = floor((length(x)-W)/O)+1;
    
    S = zeros(nfft, N_f);
    
    TimeAxis_s = ((0:N_f-1) * O + W/2) / fs;
    
    if mod(nfft, 2) == 0
        FrequencyAxis_Hz = (-nfft/2:nfft/2-1) * fs / nfft;
    else
        FrequencyAxis_Hz = (-(nfft-1)/2:(nfft-1)/2) * fs / nfft;
    end
    
    for k = 1:N_f
        StartIdx = (k - 1) * O + 1;
        StopIdx = StartIdx + W - 1;
        frame = x(StartIdx : StopIdx);
        frame = frame .* w;
        fft_frame = fftshift(fft(frame, nfft));
        S(:, k) = fft_frame;
    end
end

%% Function to get new identifier timestamp to identify recorded data uniquely
function id = getNextRunID()
    if exist('last_run_id.mat', 'file')
        load('last_run_id.mat', 'last_id');
        id = last_id + 1;
    else
        id = 1;
    end
    last_id = id;
    save('last_run_id.mat', 'last_id');
end