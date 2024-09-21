t = linspace(0,3.5,8);
% received = sin(2*pi*f*t + pi/6); received without phase diff
% received_n = sin(2*pi*f*t + pi/6 + n*pi*sin(pi/6)); with phase diff


syms m

%received = sin(2*pi*f*t + pi/6 + m*pi*sin(pi/6));
lpf_i = cos(pi/6 + pi*sin(pi/6)*m);
lpf_q = sin(pi/6 + pi*sin(pi/6)*m);

I = zeros(1,8);
Q = zeros(1,8);

for i = 0:7
    I(i+1) = subs(lpf_i, m, i);
    Q(i+1) = subs(lpf_q, m, i);
end


I_fft = fftshift(fft(I));
Q_fft = fftshift(fft(Q));

angle1 = angle(I_fft);
angle2 = angle(Q_fft);


[~, max_idx] = max(abs(I_fft));
max_phase1 = angle1(max_idx); % we expect 

% Frequency vector

f_s = 1/1; % 1 over time between samples which is 1 because m is indexed 0 to 7 in 1s
n = length(I);  % Length of the signal
frequencies = (-n/2:n/2-1)*(f_s/n);  % Frequency vector for FFT

[~, peak_index] = max((I_fft));

% Find the corresponding frequency
peak_frequency = frequencies(peak_index);

disp(peak_frequency);
theta = asin(2*pi*peak_frequency/pi);
disp(theta*180/pi);

figure;
subplot(2,1,1);
plot(t, I, 'r');
title('I (sin) in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(t, Q, 'b');
title('Q (cos) in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Plot FFT of I and Q
figure;
subplot(2,1,1);
plot(frequencies, abs(I_fft));
title('FFT of I (sin)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

subplot(2,1,2);
plot(frequencies, abs(Q_fft));
title('FFT of Q (cos)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

