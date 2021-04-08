%% --- 1. CLASSICAL AND MODERN SPECTRUM ESTIMATION --- %%

%% 1.4 Spectrum of Autoregressive Processes %%

% (a) While using the unbiased ACF estimate when finding the AR parameters,
% this estimate can be highly erratic for larger lags close to N, where
% fewer samples are available to estiamte to estimate the PSD. ACF may not
% bepositive definite, resulting in negative PSD values.
clc; clear variables; close all;

% Add data folder to path
addpath('../data');

A = [1, -2.76, 3.81, -2.65, 0.92];
x = filter(1, A, randn(1000, 1));
x = x(500:1000); leg = {'True PSD'};
[H, f] = freqz(1, A, []); P = abs(H).^2;
figure(1); hold on;
plot(f/pi, 10*log10(abs(P)));
for i = [4, 10, 14]
    [pxx, w] = pyulear(x, i); % power spectrum estimate given AR model
    plot(w/pi, 10*log10(pxx));
    leg{end+1} = sprintf('Model Order %d', i);
end
set(gca,'fontsize', 14);
xlabel('Normalized Frequency (\times \pi rad/sample)'); 
ylabel('Power/frequency (dB/(rad/sample))'); title('PSD Estimate of AR(4) Process');
legend(leg); hold off;

%% Part c
close all;
clear all;
clc;

% Add data folder to path
addpath('../data');

n_samples = [1000, 10000]; titles = {'(1000 samples)', '(10000 samples)'};
figure(1); 
for j = 1:2
    A = [1, -2.76, 3.81, -2.65, 0.92];
    x = filter(1, A, randn(n_samples(j), 1));
    x = x(500:end); leg = {'True PSD'};
    [H, f] = freqz(1, A, 256); P = abs(H).^2;
    subplot(1, 2, j); hold on;
    plot(f/pi, 10*log10(abs(P)));
    for i = [4, 10, 14]
        [pxx, w] = pyulear(x, i, 256); % power spectrum estimate given AR model
        plot(w/pi, 10*log10(pxx));
        leg{end+1} = sprintf('Model Order %d', i);
    end
    set(gca,'fontsize', 14);
    xlabel('Normalized Frequency (\times \pi rad/sample)'); 
    ylabel('Power/frequency (dB/(rad/sample))'); 
    title(strcat('PSD Estimate ',titles(j)));
    legend(leg); hold off;
end

