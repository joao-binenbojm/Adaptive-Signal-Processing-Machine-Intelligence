%% --- 1. CLASSICAL AND MODERN SPECTRUM ESTIMATION --- %%



%% 1.1. Properties of Power Spectral Density (PSD) %%
clear all;
close all;
clc;

% Simulations proof of estimate
Fs = 1000; % sampling frequency
T = 1/Fs; % sampling period
L = 3000; % length of signal
t = (0:L-1)*T; % time vector

% Pulse
S1 = zeros(size(t));
S1(1500) = 2;
Rxx1 = xcorr(S1, 'biased');
f1 = -Fs/2:Fs/length(Rxx1):Fs/2 - Fs/length(Rxx1);
P_rxx1 = abs(fftshift(fft(Rxx1)));
% Periodogram estimate
periodogram = (abs(fftshift(fft(S1))).^2)/L;
f2 = -Fs/2:Fs/length(periodogram):Fs/2 - Fs/length(periodogram);

% Sinusoidally varying signal with noise (low f)
S2 = 0.7 * sin(2*pi*0.5*t);
Rxx2 = xcorr(S2, 'biased');
P_rxx2 = abs(fftshift(fft(Rxx2)));
% Periodogram estimate
periodogram2 = (abs(fftshift(fft(S2))).^2)/L;

figure(1); 
subplot(3, 2, 1);
plot(t, S1);set(gca,'fontsize', 16); xlabel('t (s)'); ylabel('Y(t)'); title('Pulse');
subplot(3, 2, 3);
plot([-L+1:L-1], Rxx1, 'r');
set(gca,'fontsize', 14); xlabel('Time lag (k)'); ylabel('ACF'); title('Fast-decaying ACF')
subplot(3, 2, 5); hold on;
plot(f1, P_rxx1, 'r'); plot(f2, periodogram, '--b');
set(gca,'fontsize', 14); xlabel('Frequency (Hz)'); ylabel('Magnitude'); title('PSD');
legend('Def. 1', 'Def. 2'); ylim([0 0.002]); hold off;
%%%%%%
subplot(3, 2, 2);
plot(t,S2); set(gca,'fontsize', 14); xlabel('t (s)'); ylabel('Y(t)'); title('Sinusoidal signal');
subplot(3, 2, 4);
plot([-L+1:L-1], Rxx2);
set(gca,'fontsize', 14); 
xlabel('Time lag (k)'); ylabel('ACF'); title('Slow-decaying ACF');
subplot(3, 2, 6); hold on;
plot(f1, P_rxx2, 'r');
plot(f2, periodogram2, '--b');
set(gca,'fontsize', 14); 
xlabel('Frequency (Hz)'); ylabel('Magnitude'); title('PSD')
legend('Def. 1', 'Def. 2');
hold off;


