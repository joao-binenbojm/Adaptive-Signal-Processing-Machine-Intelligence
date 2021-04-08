%% --- 1. CLASSICAL AND MODERN SPECTRUM ESTIMATION --- %%


%% 1.2 Periodogram-based methods Applied to Real-World Data %%
close all;
clear all;
clc;

% Add data folder to path
addpath('../data');

% a. Sunspot time-series
load sunspot.dat

figure(1); plot(sunspot(:, 1), sunspot(:, 2));
xlabel('Time (years)'); ylabel('Sunspot Count'); title('Sunspot Time Series');
data = sunspot(:, 2);

% Non-log data
[P_raw, w] = periodogram(data); % PSD for raw data
P_raw(P_raw < 0.001) = 1;
w = w/pi;
P_centered = periodogram(data - mean(data)); % PSD for centered data
P_centered(P_centered < 0.01) = 1;
P_detrend = periodogram(detrend(data)); % PSD for detrended, centered data
P_detrend(P_detrend < 0.01) = 1;

% Log data
data = log10(data + eps);
P_raw_log = periodogram(data); % PSD for raw data
P_raw_log(P_raw_log < 0.01) = 1;
P_centered_log = periodogram(data - mean(data)); % PSD for centered data
P_centered_log(P_centered_log < 0.01) = 1;
P_detrend_log = periodogram(detrend(data)); % PSD for detrended, centered data
P_detrend_log(P_detrend_log < 0.01) = 1;

% Plotting PSD estimates
figure(2);
subplot(1, 2, 1); hold on;
plot(w, 10*log10(P_raw + eps)); 
plot(w, 10*log10(P_centered + eps));
plot(w, 10*log10(P_detrend + eps));
set(gca,'fontsize', 14);
xlabel('Normalized Frequency (\times \pi rad/sample)'); 
ylabel('Power/frequency (dB/(rad/sample))'); title('PSD Estimates (Original Data)');
legend('Raw', 'Centered', 'Detrended'); hold off;
% Log transformed data
subplot(1, 2, 2); hold on;
plot(w, 10*log10(P_raw_log + eps)); 
plot(w, 10*log10(P_centered_log + eps));
plot(w, 10*log10(P_detrend_log + eps));
set(gca,'fontsize', 14);
xlabel('Normalized Frequency (\times \pi rad/sample)'); 
ylabel('Power/frequency (dB/(rad/sample))'); title('PSD Estimates (Log Data)');
legend('Raw', 'Centered', 'Detrended'); hold off;
%% b. EEG
clc; clear variables; close all;

% Add data folder to path
addpath('../data');

load EEG_Data_Assignment1.mat;
f_stimulus_range = [11:20]; T = 1/fs; L = length(POz);
POz = normalize(POz);
[P_c, f] = pwelch(POz, L, 0, L, fs); % Ensures 5 DFT samples/Hz
wl = [1, 5, 10]/T; % Number of samples per window

% Plotting PSD Estimates
figure(1); subplot(2, 2, 1);
plot(f(1:round(length(f)/6)), 10*log10(P_c(1:round(length(f)/6))));
set(gca,'fontsize', 14);
xlabel('Frequency (Hz)');
ylabel('Power/frequency (dB/(rad/sample))'); title('EEG PSD');
for i=1:3
    [P_mean, f] = pwelch(POz, wl(i), 0, L, fs);
    subplot(2, 2, i+1);
    plot(f(1:round(length(f)/6)), 10*log10(P_mean(1:round(length(f)/6)))); 
    set(gca,'fontsize', 14); xlabel('Frequency (Hz)');
    ylabel('Power/frequency (dB/(rad/sample))'); title(sprintf('EEG PSD (W=%d s)', wl(i)*T));
end

% Pretty bad move to make the window length too small, since it makes the
% estimate have less variance, yet more biased. From this example, we can't
% even see the SSVEP anymore.
