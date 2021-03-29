%% --- 1. CLASSICAL AND MODERN SPECTRUM ESTIMATION --- %%

%% 1.5 Real World Signals: Respiratory Sinus Arrhythmia from RR-Intervals %%

%% Part a
close all; clear all; clc;
load ECG_data;

xRRI_1 = detrend(normalize(xRRI_1));
xRRI_2 = detrend(normalize(xRRI_2));
xRRI_3 = detrend(normalize(xRRI_3));

RRI_data = {xRRI_1; xRRI_2; xRRI_3};
leg = {'Standard', '(W = 50s)', '(W = 150s)'};
figure(1);
for j = 1:3 % Trials
    subplot(1, 3, j); hold on;
    L = length(RRI_data{j}); wl = [L, 50*fsRRI_1, 150*fsRRI_1];
    for i = 1:3 % Window sizes
        [P_mean, w] = pwelch(RRI_data{j}, wl(i), 0, L, 4);
        plot(w, 10*log10(P_mean));
    end
    set(gca,'fontsize', 14);
    xlabel('Frequency (Hz)');
    ylabel('Power/frequency (dB/(rad/sample))'); 
    title(sprintf('RRI PSD Estimate (Trial %d)', j));
    legend(leg); hold off;
end

% (b) The first trial for unconstrained breathing gives us no discernable
% frequencies, as the subject was allowed to breathe at any rhythm. There
% were discernable peaks in both trial 2 and 3 at ?

%% Part c
clc; clear all; close all;
load ECG_data;

xRRI_1 = detrend(normalize(xRRI_1));
xRRI_2 = detrend(normalize(xRRI_2));
xRRI_3 = detrend(normalize(xRRI_3));

% 'Determining model order' code
% leg = {};
% figure(1); hold on;
% for i = [1:4]
%     [pxx, w] = pyulear(xRRI_1, i, 2048); % power spectrum estimate given AR model
%     plot(w/pi, 10*log10(pxx));
%     leg{end+1} = sprintf('Model Order %d', i);
% end
% xlabel('Normalized Frequency (\times \pi rad/sample)');
% ylabel('Power/frequency (dB/(rad/sample))'); 
% legend(leg); hold off;

% Trial 1 -> p=3, Trial 2 -> p=7, Trial 3 -> p=4

% Plotting results
order = [3, 7, 4]; RRI_data = {xRRI_1; xRRI_2; xRRI_3};
figure(1); hold on;
for i = [1, 2, 3]
    [pxx, w] = pyulear(RRI_data{i}, order(i), 2048, 4); % power spectrum estimate given AR model
    plot(w, 10*log10(pxx));
end
set(gca,'fontsize', 14);
xlabel('Frequency (Hz)');
ylabel('Power/frequency (dB/(rad/sample))'); 
title('AR PSD Estimate of RRI Data');
legend('Trial 1', 'Trial 2', 'Trial 3'); hold off;

%% Confirming frequencies with MUSIC algorithm
clc; clear all; close all;
load ECG_data;

xRRI_1 = detrend(normalize(xRRI_1));
xRRI_2 = detrend(normalize(xRRI_2));
xRRI_3 = detrend(normalize(xRRI_3));
RRI_data = {xRRI_1; xRRI_2; xRRI_3}; order = [3, 7, 4];

for i = [1, 2, 3]
    [X,R] = corrmtx(RRI_data{i}, 14, 'mod'); % each row of X is a separate row of the correlation matrix
    [S,F] = pmusic(R, 4, [ ], 4, 'corr');
    figure(i);
    plot(F, S); xlabel('Frequency (Hz)');
    ylabel('Pseudospectrum'); title (sprintf('Trial %d', i));
    
end

