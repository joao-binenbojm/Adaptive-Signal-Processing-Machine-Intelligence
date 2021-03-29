%% --- 2. Adaptive signal processing --- %%

%% 2.2. Adaptive Noise Cancellation %%

%% (a) Implementing ALE
clc; clear variables; close all;
M = [5:5:20]; delta = [3:25]; lr = 0.01;
eta = filter([1, 0, 0.5], [1], randn(1500, 1));
eta = eta(501:end);
x = sin([1:1000]*0.01*pi)';
s = x + eta;

MSPE = zeros(length(M), length(delta));
L = length(s); leg = {};

figure(1);
% mspe vs order&delay
subplot(1, 2, 1); hold on; 
set(gca,'fontsize', 16); 
for i = 1:length(M)
    for j = 1:length(delta)
        [~,xhat,~] = ale_lms(s, lr, delta(j), M(i));
        MSPE(i, j) = (1/L)*(x-xhat)' * (x-xhat);
    end
    plot(delta, 10*log10(MSPE(i,:)));
    leg{i} = ['M = ' num2str(M(i))];
end
legend(leg); xlim([3, 25]);
xlabel('Delay ($\Delta$)', 'Interpreter', 'Latex');
ylabel('MSPE (dB)'); title('MSPE vs. Filter Order & Delay');

% Showing the improvement in obtaining noisy signal
[w, xhat, error] = ale_lms(s, 0.01, 3, 5);
subplot(1, 2, 2); hold on;
set(gca,'fontsize', 16); 
plot([1:length(s)], s, 'k');
plot([1:length(xhat)], xhat, 'g');
plot([1:length(x)], x, 'r', 'LineWidth', 2); 
legend('Noisy signal', 'ALE', 'Clean signal')
xlabel('Sample (n)'); ylabel('Signal'); title('ALE Signal Denoising');
hold off;



%% (b) Testing out ALE system
clc; clear variables; close all;
lr = 0.01;
eta = filter([1, 0, 0.5], [1], randn(1500, 1));
eta = eta(501:end);
x = sin([1:1000]*0.01*pi)';
s = x + eta;
L = length(s);

M = 5; delta = 3;
MSPE = zeros(2, 2);

figure(1);
% mspe vs order&delay
for i = 1:20
    % Filter order MSPE
    [~,xhat,~] = ale_lms(s, lr, delta, i);
    MSPE(1, i) = (1/L)*(x-xhat)' * (x-xhat);
    % Delay plot MSPE
    [~,xhat,~] = ale_lms(s, lr, i, M);
    MSPE(2, i) = (1/L)*(x-xhat)' * (x-xhat);
end
% Filter order
subplot(1, 2, 1); hold on; set(gca,'fontsize', 16);
plot([1:20], 10*log10(MSPE(1,:)), 'r');  hold off;
xlabel('Filter Order (M)'); ylabel('MSPE (dB)');
title('MSPE vs. Filter Order ($\Delta=3$)', 'Interpreter', 'Latex'); 
% Delay
subplot(1, 2, 2); hold on; set(gca,'fontsize', 16); 
plot([1:20], 10*log10(MSPE(2,:))); hold off;
xlabel('Delay ($\Delta$)', 'Interpreter', 'Latex'); 
ylabel('MSPE (dB)');
title('MSPE vs. Delay (M=5)', 'Interpreter', 'Latex'); 


%% (c) Testing out ANC system
clc; clear variables; close all;
lr = 0.01;
eta = filter([1, 0, 0.5], [1], randn(1500, 1));
eta = eta(501:end);
x = sin([1:1000]*0.01*pi)';
s = x + eta;
sec_noise = 0.5*eta - 0.1*delayseq(eta, 2); % delay of 2
L = length(s);

MSPE = zeros(2, 2);

[~,xhat] = ale_lms(s, lr, 3, 5);
[~,xhat_anc] = anc_lms(s, sec_noise, lr, 10);

% Signal reconstruction via noise cancellation
figure(1); subplot(1, 2, 1);
hold on; set(gca,'fontsize', 16);
plot([1:length(xhat)], xhat, 'k');
plot([1:length(xhat_anc)], xhat_anc, 'g');
plot([1:length(x)], x, 'red');
xlabel('Sample (n)'); ylabel('Signal'); 
title('ALE vs ANC For Signal Reconstruction');
legend('ALE', 'ANC', 'Clean');
hold off;

% Plotting Computed MSPE dynamics
MSPE_ale = 10*log10(movmean((x-xhat).^2, 30));
MSPE_anc = 10*log10(movmean((x-xhat_anc).^2, 30));

subplot(1, 2, 2);
hold on; set(gca,'fontsize', 16);
plot([1:length(MSPE_ale)], MSPE_ale, 'r');
plot([1:length(MSPE_anc)], MSPE_anc, 'b');
xlabel('Sample (n)'); ylabel('MSPE (dB)'); 
title('MSPE: ALE vs ANC');
legend('ALE', 'ANC');
hold off;

%% (d) EEG Data!! 
clc; clear variables; close all;
load EEG_Data_Assignment2;

% Normalizing our data
% Cz_norm = normalize(Cz);
% data = detrend(POz);
data = detrend(Cz);

% figure(2); subplot(1, 2, 1);
% ylim([0, 55]); xlim([0, 4.8]); hold on; 
% set(gca,'fontsize', 16);
% spectrogram(data, rectwin(3540), round(0.3*(3540)), 16382, fs, 'yaxis');
% title('Spectrogram of Original Cz Data', 'Interpreter', 'Latex');
% hold off;

t = [0:length(data)-1];
lr = [0.01, 0.001, 0.0001];
M = [1, 5, 10, 20];
stdevs = [2];
for i = stdevs
    mains = sin((2*pi*50/fs)*t) + (10^(-i))*randn(1, length(data));
    mains = mains';
    counter = 1;
    figure(i);
    sprintf('Variance: 1e-%d',i)
    for j = 1:length(lr)
        for k = 1:length(M)
            [w,xhat] = anc_lms(Cz, mains, lr(j), M(k));
%             xhat = anc_denoise(POz, mains, mean(w(:,end-100:end),2));
%             subplot(1, 2, 2);
            subplot(length(lr), length(M), counter)
            ylim([0, 55]); xlim([1, 4.8]); hold on; set(gca,'fontsize', 16);
            spectrogram(xhat, rectwin(3540), round(0.3*(3540)), 16382, fs, 'yaxis');
            title(sprintf('i=%d, lr=%.4f, M=%d ', i, lr(j), M(k)));
%             title('Spectrogram of Cleaned Cz Data (M=5, $\mu=0.01$)', 'Interpreter', 'Latex');
            hold off;
            counter = counter + 1;
        end
    end
end