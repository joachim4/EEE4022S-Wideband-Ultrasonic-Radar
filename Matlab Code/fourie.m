% Audio Radar

%% Clear variables and command window
clear all;
close all;
clc;

%% Define constants and parameters
SpeedSoundWave_ms = 343;             % [m/s]  ->speed of wave
Fc_Hz = 40000;                        % [Hz]
TimeDuration_s = 10;                  % [s]  
Fs = 192000;                         % [Hz]   ->Sampling rate. So 44 100 samples are obtained in one second

%% 
Ts = 1/Fs;                           % Sampling periodf
t = 0:Ts:(TimeDuration_s);           % time vector for pulse

%% Generate the transmit signal
% Pure sinusoid 
TxSignal = 1*sin(2*pi*Fc_Hz*t);

%% Play out transmit signal through the speakers
soundsc(TxSignal, Fs, 24) % Transmit the signal

%% Record received samples from the microphone
RecLength_samples = length(TxSignal);
RecLength_s = RecLength_samples*1/Fs; 
recObj = audiorecorder(Fs, 24, 1);
recordblocking(recObj, RecLength_s);  % Records audio for a fixed number of seconds
RX_signal = getaudiodata(recObj);   % Store recorded audio signal in double-precision array

%% Plot the time-domain signals
figure; axes('fontsize', 12);
subplot(2,1,1);
plot(t, TxSignal); % plot transmit signal
xlabel('Time (s)', 'fontsize', 12);
ylabel('Amplitude (linear)', 'fontsize', 12);
title('Transmit signal', 'fontsize', 12);
grid on;

subplot(2,1,2);
plot(t, RX_signal); % plot received signal
xlabel('Time (s)', 'fontsize', 12);
ylabel('Amplitude (linear)', 'fontsize', 12);
title('Received signal', 'fontsize', 12);
grid on;

%% FFT and FFTShift

% Compute the FFT of the signals
N = length(t);
f = (-N/2:N/2-1)*(Fs/N);  % Frequency vector for plotting

TxSignal_FFT = fft(TxSignal);
TX_signal_FFTShifted = fftshift(TxSignal_FFT);

RX_signal_FFT = fft(RX_signal);
RX_signal_FFTShifted = fftshift(RX_signal_FFT);

% Plot the FFT shifted signals
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
