%% This code was created by Joachim Gengan GNGJOA003 - in partial fulfulment of EEE4022S research project
%% The code has been adapted from example code provided by Dr M. Y. Abdul Gaffar

close all;

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

lamda = SpeedSoundWave_ms/Fc_Hz;

%% Unique identifier for this run using current timestamp on computer to differentiate plots
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
run_id = getNextRunID();

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

%% Saving figures for later use to be plot for report
%saveas(gcf, ['Figures\TransmitReceiveSignals\TransmitReceiveSignals_' timestamp '_' num2str(run_id) '.png']);
%savefig(['Figures\TransmitReceiveSignals\TransmitReceiveSignals_' timestamp '_' num2str(run_id) '.fig']);

%% Digital filtering - removes unwanted frequencies outside of 20 kHz - 60 kHz passband if required, also can apply a notch (sharp stop band filter at 40 kHz) 

lowFreq = 20000;  
highFreq = 60000;
[b_bp, a_bp] = butter(6, [lowFreq highFreq]/(Fs/2), 'bandpass');
RX_signal_bp = filtfilt(b, a, RX_signal);

w0 = Fc_Hz/(Fs/2);  
bw = 0.003;  % Bandwidth of the notch filter at 40 kHz   
[b_notch, a_notch] = iirnotch(w0, bw); 
RX_signal_filtered = filtfilt(b_notch, a_notch, RX_signal);

%% Frequency response of the notch filter
figure;
freqz(b_notch, a_notch, 1024, Fs); % General frequency response
title('Frequency Response of the Notch Filter');
grid on;

%% Notch filter response
figure;
[H, f] = freqz(b_notch, a_notch, 1024, Fs);
plot(f, 20*log10(abs(H)));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Zoomed Frequency Response of the Notch Filter');
grid on;

% View the frequency response of the notch filter
xlim([Fc_Hz-500 Fc_Hz+500]);  % Frequency limits for better visualisation
ylim([-50 5]); 


%% Downmixing - brings signal at high frequency down to baseband
cos_component = cos(2 * pi * Fc_Hz * t);
sin_component = -sin(2 * pi * Fc_Hz * t);
I_signal = RX_signal_filtered .* cos_component;
Q_signal = RX_signal_filtered .* sin_component;

fcutoff = 20000;
[b, a] = butter(5, fcutoff / (Fs / 2), 'low');

I_baseband = filter(b, a, I_signal);
Q_baseband = filter(b, a, Q_signal);
baseband_signal = I_baseband + (1i * Q_baseband);

%% STFT calculation
% Parameters for STFT 
W = 1024; % Frame length 
O = 0.5 * W; %  overlap
nfft = 1024; % FFT length
w = hamming(W);

[S, TimeAxis_s, FrequencyAxis_Hz] = JoachimSpectrogram(baseband_signal, w, W, O, Fs, nfft);

maxFreq = 6000; % Can set limits for the Doppler frequency graph so that it will be more zoomed in
minFreq = -6000;
freq_idx = find((FrequencyAxis_Hz <= maxFreq) & (FrequencyAxis_Hz >= minFreq));
FreqVectorOfInterest = FrequencyAxis_Hz(freq_idx);
F_OfInterest = S(freq_idx, :);
F_OfInterestToPlot = abs(F_OfInterest)/max(max(abs(F_OfInterest))); % normalise plot

% Creates a spectrogram plot here
clims = [-40 0]; % Sets the dynamic range of the spectrogram
figure;
%imagesc(TimeAxis_s, FrequencyAxis_Hz/1e3, 20*log10(abs(S)),clims);
imagesc(TimeAxis_s,FreqVectorOfInterest,20*log10(F_OfInterestToPlot), clims);
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram of Received Signal');
colormap('jet');
axis xy;

%saveas(gcf, ['Figures3\Spectrograms\Spectrogram_' timestamp '_' num2str(run_id) '.png']);
%savefig(['Figures3\Spectrograms\Spectrogram_' timestamp '_' num2str(run_id) '.fig']);

maxSpeed_m_s = 25; % Can set limits for the m/s graph so that it will be more zoomed in
minSpeed_m_s = -25; 

maxSpeed_km_hr = 70;
minSpeed_km_per_hr = -70; 
speed_m_per_sec = FrequencyAxis_Hz*lamda/2;
speed_km_per_hr = speed_m_per_sec*(60*60/1000);
speed_km_per_hr_Idx = find((speed_km_per_hr <= maxSpeed_km_hr) & (speed_km_per_hr >= minSpeed_km_per_hr));
SpeedVectorOfInterest = speed_km_per_hr(speed_km_per_hr_Idx);

S_OfInterest = S(speed_km_per_hr_Idx, :);
S_OfInterestToPlot = abs(S_OfInterest)/max(max(abs(S_OfInterest))); % normalise plot

%speed_m_per_sec = FrequencyAxis_Hz*lamda/2;
speed_m_per_sec_Idx = find((speed_m_per_sec <= maxSpeed_m_s) & (speed_m_per_sec >= minSpeed_m_s));
SpeedVectorOfInterest_m_s = speed_m_per_sec(speed_m_per_sec_Idx);

S_OfInterest_m_s = S(speed_m_per_sec_Idx, :);
S_OfInterestToPlot_m_s = abs(S_OfInterest_m_s)/max(max(abs(S_OfInterest_m_s))); % normalise plot

% Plot of spectrogram in kilometres per hour
clims = [-35 0];
figure; imagesc(t,SpeedVectorOfInterest,20*log10(S_OfInterestToPlot), clims);
xlabel('Time (s)');
ylabel('Speed (km/hr)');
title('Spectrogram: zoomed in (km/hr)');
grid on;
colorbar;
colormap('jet');
axis xy;

%saveas(gcf, ['Figures3\Spectrograms\SpectrogramKMPH_' timestamp '_' num2str(run_id) '.png']);
%savefig(['Figures3\Spectrograms\SpectrogramKMPH_' timestamp '_' num2str(run_id) '.fig']);

% Plot of spectrogram in meters per second
clims = [-35 0];
figure; imagesc(t,SpeedVectorOfInterest_m_s,20*log10(S_OfInterestToPlot_m_s), clims);
xlabel('Time (s)');
ylabel('Speed (m/s)');
%title('Spectrogram: zoomed in (m/s)');
grid on;
colorbar;
colormap('jet');
axis xy;

% Normalisation 2 - helps to see instantaneous speed better, most useful for targets with micro-Doppler signatures like a human walking
S_OfInterestToPlot_m_s_n2 = abs(S_OfInterest_m_s)./max(abs(S_OfInterest_m_s)); 

% Plot of spectrogram in meters per second
clims = [-15 0];
figure; imagesc(t,SpeedVectorOfInterest_m_s,20*log10(S_OfInterestToPlot_m_s_n2), clims);
xlabel('Time (s)');
ylabel('Speed (m/s)');
title('Spectrogram: zoomed in (m/s)');
grid on;
colorbar;
colormap('jet');
axis xy;

%saveas(gcf, ['Figures3\Spectrograms\SpectrogramMS_' timestamp '_' num2str(run_id) '.png']);
%savefig(['Figures3\Spectrograms\SpectrogramKMS_' timestamp '_' num2str(run_id) '.fig']);

%save(['Workspaces3\SpectrogramWorkspace_' timestamp '_' num2str(run_id) '.mat']); 

%% Custom STFT function
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

%% Function to get new identifier to store data as a unique entry for easily finding plots later
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