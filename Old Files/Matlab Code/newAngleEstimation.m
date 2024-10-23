% Define time vector
t = linspace(0, 3.5, 8);

% Define symbolic variable for phase difference
syms m

% Define I and Q components
lpf_i = cos(pi/6 + pi*sin(pi/6)*m);
lpf_q = cos(pi/6 + pi*sin(pi/6)*m);

% Initialize I and Q arrays
I = zeros(1,8);
Q = zeros(1,8);

% Calculate I and Q values
for i = 0:7
    I(i+1) = subs(lpf_i, m, i);
    Q(i+1) = subs(lpf_q, m, i);
end

% Perform FFT on I and Q
I_fft = fftshift(fft(I));
Q_fft = fftshift(fft(Q));

% Calculate angles
angle1 = angle(I_fft);
angle2 = angle(Q_fft);

% Find peak frequency
[~, max_idx] = max(abs(I_fft));
max_phase1 = angle1(max_idx);

% Calculate frequency vector
f_s = 1; % Sampling frequency (1 over time between samples)
n = length(I);
frequencies = (-n/2:n/2-1)*(f_s/n);

% Find peak frequency
[~, peak_index] = max(abs(I_fft));
peak_frequency = frequencies(peak_index);

% Calculate theta
theta = asin(2*pi*peak_frequency/pi);

% Display results
disp(['Peak Frequency: ', num2str(peak_frequency)]);
disp(['Theta (degrees): ', num2str(theta*180/pi)]);

% Plot time domain signals
figure;
subplot(2,1,1);
plot(t, I, 'r');
title('I (cos) in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(t, Q, 'b');
title('Q (sin) in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Plot FFT of I and Q
figure;
subplot(2,1,1);
plot(frequencies, abs(I_fft));
title('FFT of I (cos)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

subplot(2,1,2);
plot(frequencies, abs(Q_fft));
title('FFT of Q (sin)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Perform cross-correlation between I and Q
[correlation, lags] = xcorr(I, Q);

% Find the peak of the correlation
[~, max_corr_idx] = max(abs(correlation));
lag = lags(max_corr_idx);

% Determine which signal leads
if lag > 0
    disp('I signal leads Q signal');
elseif lag < 0
    disp('Q signal leads I signal');
else
    disp('I and Q signals are in phase');
end

% Plot cross-correlation
figure;
plot(lags, correlation);
title('Cross-correlation between I and Q');
xlabel('Lag');
ylabel('Correlation');
grid on;

% Calculate and display the phase difference
phase_diff = atan2(Q(1), I(1)) * 180 / pi;
disp(['Initial phase difference (Q relative to I): ', num2str(phase_diff), ' degrees']);

% Calculate the phase difference directly
phase_diff_direct = unwrap(angle(Q + 1i*I)) * 180 / pi;

% Calculate the average phase difference
avg_phase_diff = mean(diff(phase_diff_direct));
disp(['Average phase difference: ', num2str(avg_phase_diff), ' degrees']);

% Plot the phase difference
figure;
plot(t, phase_diff_direct);
title('Phase Difference (Q relative to I)');
xlabel('Time (s)');
ylabel('Phase Difference (degrees)');
grid on;

% Plot I vs Q (phasor diagram)
figure;
plot(I, Q, 'b-', 'LineWidth', 2);
hold on;
plot(I, Q, 'ro', 'MarkerSize', 8);
title('I vs Q (Phasor Diagram)');
xlabel('I');
ylabel('Q');
axis equal;
grid on;