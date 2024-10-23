%% Clear variables and command window
clear all;
close all;
clc;

%% Define constants and parameters
SpeedSoundWave_ms = 343; % [m/s] -> Speed of sound wave
Fc_Hz = 40000; % [Hz] -> Carrier frequency
TimeDuration_s = 10; % [s] -> Duration of signal transmission and reception
Fs = 192000; % [Hz] -> Sampling rate
v_target_kmh = 50; % [km/h] -> Velocity of the target
v_target_ms = v_target_kmh / 3.6; % [m/s] -> Convert velocity to m/s

%% Calculate Doppler Frequency
f_doppler = (2 * v_target_ms * Fc_Hz) / SpeedSoundWave_ms;
disp(['Calculated Doppler frequency: ' num2str(f_doppler) ' Hz']);

%% Generate the transmit signal
Ts = 1/Fs; % Sampling period
t = 0:Ts:TimeDuration_s;
TxSignal = sin(2 * pi * Fc_Hz * t);

%% Play out transmit signal through the speakers
soundsc(TxSignal, Fs, 24); % Transmit the signal

%% Record received samples from the microphone
RecLength_samples = length(TxSignal);
RecLength_s = RecLength_samples * Ts;
recObj = audiorecorder(Fs, 24, 1);
recordblocking(recObj, RecLength_s); % Record audio for a fixed number of seconds
RX_signal = getaudiodata(recObj); % Store recorded audio signal in double-precision array

%% Save the recorded signal to a file
save('received_signal.mat', 'RX_signal', 'Fs', 'Fc_Hz', 't');