%% --- 1. CLASSICAL AND MODERN SPECTRUM ESTIMATION --- %%

%% 1.3 Correlation Estimation %%
close all;
clear all;
clc;

% Add data folder to path
addpath('../data');

% Code validation
t = [0:0.01:9.99]; fs = 100; x1 = 6*sin(2*pi*20*t) + 2*randn(size(t)); x2 = randn(size(t));
model = arima('Constant', 0,'AR',{2.21, -2.94, 2.17, -0.96},'Variance',1); x3 = simulate(model, length(t));
[P1, f] = correlogram(x1, 'biased', fs); [P2, f] = correlogram(x2, 'biased', fs); [P3, f] = correlogram(x3, 'biased', fs);
[P1_unbiased, f] = correlogram(x1, 'unbiased', fs); [P2_unbiased, f] = correlogram(x2, 'unbiased', fs); [P3_unbiased, f] = correlogram(x3, 'unbiased', fs);
[real_p, f_new] = periodogram(x1, [], length(x1), fs);

f = 2*f/fs; % Going from Hz to normalized frequency

% % ACFs
figure(1);
% Biased
subplot(2, 3, 1); plot(t, acf(x1, 'biased')); 
set(gca,'fontsize', 14); xlabel('Time (s)'); 
ylabel('ACF'); title('ACF of Noisy Sinusoid (Biased)');
subplot(2, 3, 2); plot(t, acf(x2, 'biased')); 
set(gca,'fontsize', 14); xlabel('Time (s)');
ylabel('ACF'); title('ACF of WGN (Biased)');
subplot(2, 3, 3); plot(t, acf(x3, 'biased')); 
set(gca,'fontsize', 14); xlabel('Time (s)');
ylabel('ACF'); title('ACF of AR(4) Process (Biased)');
% Unbiased
subplot(2, 3, 4); plot(t, acf(x1, 'unbiased')); 
set(gca,'fontsize', 14); xlabel('Time (s)'); 
ylabel('ACF'); title('ACF of Noisy Sinusoid (Unbiased)');
subplot(2, 3, 5); plot(t, acf(x2, 'unbiased')); 
set(gca,'fontsize', 14); xlabel('Time (s)');
ylabel('ACF'); title('ACF of WGN (Unbiased)');
subplot(2, 3, 6); plot(t, acf(x3, 'unbiased')); 
set(gca,'fontsize', 14); xlabel('Time (s)');
ylabel('ACF'); title('ACF of AR(4) Process (Unbiased)');
%PSDs
figure(2);
% Biased
subplot(3, 2, 1); plot(f, 10*log10(P1)); 
set(gca,'fontsize', 14); xlabel('Normalized Frequency (\times \pi rad/sample)'); 
ylabel('Pow/freq (dB/(rad/sample))'); title('PSD of Noisy Sinusoid (Biased)');
subplot(3, 2, 3); plot(f, 10*log10(P2)); 
set(gca,'fontsize', 14); xlabel('Normalized Frequency (\times \pi rad/sample)');
ylabel('Pow/freq (dB/(rad/sample))'); title('PSD of WGN (Biased)');
subplot(3, 2, 5); plot(f, 10*log10(P3)); 
set(gca,'fontsize', 14); xlabel('Normalized Frequency (\times \pi rad/sample)');
ylabel('Pow/freq (dB/(rad/sample))'); title('PSD of AR(4) Process (Biased)');
% Unbiased
subplot(3, 2, 2); plot(f, 10*log10(P1_unbiased)); 
set(gca,'fontsize', 14); xlabel('Normalized Frequency (\times \pi rad/sample)'); 
ylabel('Pow/freq (dB/(rad/sample))'); title('PSD of Noisy Sinusoid (Unbiased)');
subplot(3, 2, 4); plot(f, 10*log10(P2_unbiased)); 
set(gca,'fontsize', 14); xlabel('Normalized Frequency (\times \pi rad/sample)');
ylabel('Pow/freq (dB/(rad/sample))'); title('PSD of WGN (Unbiased)');
subplot(3, 2, 6); plot(f, 10*log10(P3_unbiased)); 
set(gca,'fontsize', 14); xlabel('Normalized Frequency (\times \pi rad/sample)');
ylabel('Pow/freq (dB/(rad/sample))'); title('PSD of AR(4) Process (Unbiased)');

%% Part b/c 
clc; clear all; close all;

% Add data folder to path
addpath('../data');

t = [0:0.01:9.99]; fs = 100; P1_store = [];
noisy_sinusoid = repmat(sin(2*pi*20*t) - sin(2*pi*10*t), 1000, 1)' + 10*randn(1000, 1000);
figure(1);
for i = 1:100
    [P1, f] = correlogram(noisy_sinusoid(:, i), 'biased', fs);
    f = 2*f/fs;
    subplot(2, 2, 1); hold on;
    plot(f, P1, 'b'); hold off;
    subplot(2, 2, 3); hold on;
    plot(f, 10*log10(P1), 'b'); hold off;
    P1_store = [P1_store; P1];  %P2_store = [P2_store; P2];
end
% WGN
subplot(2, 2, 1); hold on; plot(f, mean(P1_store, 1), 'y'); set(gca,'fontsize', 14);
xlabel('Normalized Frequency (\times \pi rad/sample)'); ylabel('Power/frequency (rad/sample)');
title('PSD estimate noisy sinusoid'); legend('Realisations', 'Mean'); hold off;
subplot(2, 2, 2); hold on; plot(f, std(P1_store, 1), 'r'); set(gca,'fontsize', 14);
xlabel('Normalized Frequency (\times \pi rad/sample)'); ylabel('Power/frequency (rad/sample)');
title('Standard deviation of PSD estimate'); hold off;
% Noisy sinusoid
subplot(2, 2, 3); hold on; plot(f, 10*log10(mean(P1_store, 1)), 'y'); set(gca,'fontsize', 14);
xlabel('Normalized Frequency (\times \pi rad/sample)'); ylabel('Power/frequency (dB/(rad/sample))');
title('PSD estimate of noisy sinusoid (dB)'); legend('Realisations', 'Mean');hold off;
subplot(2, 2, 4); hold on; plot(f, 10*log10(std(P1_store, 1)), 'r'); set(gca,'fontsize', 14);
xlabel('Normalized Frequency (\times \pi rad/sample)'); ylabel('Power/frequency (dB/(rad/sample))');
title('Standard deviation of PSD estimate (dB)'); hold off;



%% Part d - Generating complex exponentials
clc; clear all; close all;

% Add data folder to path
addpath('../data');

t = [0:0.02:9.98]; fs = 50; 
x = 2*exp(i * 2 * pi * 19 * t) + 3*exp(i * 2 * pi * 20 * t);
[P1, f1] = correlogram(x(1:20), 'biased', fs);
[P2, f2] = correlogram(x(1:80), 'biased', fs);
[P3, f3] = correlogram(x, 'biased', fs);
f1 = 2*f1/fs; f2 = 2*f2/fs; f3 = 2*f3/fs;

figure(1); hold on;
plot(f1, 10*log10(P1), 'r'); plot(f2, 10*log10(P2), 'b'); plot(f3, 10*log10(P3), 'g');
set(gca,'fontsize', 14); xlabel('Normalized Frequency (\times \pi rad/sample)'); 
ylabel('Power/frequency (dB/(rad/sample))');
legend('20 samples', '80 samples', '500 samples'); 
title('PSD of complex exponential');  hold off; xlim([0.6 1]); 

%% Part e - MUSIC method
clc; clear all; close all;
n = 0:30;
noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+ noise;
tic
[X,R] = corrmtx(x, 14, 'mod'); % each row of X is a separate row of the correlation matrix
[S,F] = pmusic(R, 2, [ ], 1, 'corr'); % p=dim(signal_subspace), 'corr' means that we interpret X as autocorrealation matrix.
time_elapsed1 = toc
tic
[Pxx, w] = periodogram(x);
time_elapse2 = toc
figure(1); subplot(1, 2, 1);
plot(2*F,S,'linewidth', 2); set(gca,'xlim',[0.50 0.80]); set(gca,'fontsize', 14);
grid on; xlabel('Normalized Frequency (\times \pi rad/sample)'); 
ylabel('Pseudospectrum'); title('MUSIC Pseudospectrum')
subplot(1, 2, 2); plot(w/pi, Pxx); set(gca,'xlim',[0.50 0.80]); set(gca,'fontsize', 14);
grid on; ylabel('Power/frequency (dB/(rad/sample))');
xlabel('Normalized Frequency (\times \pi rad/sample)'); title('Standard Periodogram Estimate')


% corrmtx - Give a length N vector S and a positive integer M,
% return an N+M by M+1 matrix X such that X' * X is a (biased) estimate
% of the M+1 by M+1 autocorrelation matrix of S (N > M).
% M is the prediction model order

% Method is how autocorrelation is computed.'mod' means that X is (N-M)by(M+1)
% rectangular Toeplitz matrix tha generate an autocorrelation estimate for
% the length-n data vector x, derived using forward and backward prediction
% error estimates. Matrix can be used to perform autoregressive parameter
% estimation using th emoedified covariance method.

% R is what we want, the M+1 by M+1 correlation matrix R (~ X'*X). This
% method fits an AR model to the signal!! uses forward-backward method.
% AR parameters are used for autocorrelation estimate.

% pmusic --> Implements the Multiple Signal Classification algorithm for
% frequency estimation. substantial performance advantages, yet achieved
% for a given cost in computation (searching in parameter space) and
% storage (of array calibration data). Uses eigenspace method.
% Since Rxx is Hermitian, all of its M eigenvectors are orthogonal. 
% vectors corresponding to p largest eigenvalues span the signal subspace
% Us.  The remaining ones span the noise subspace, Un
% MUSIC outperforms simple methods such as picking peaks of DFT spectra in
% the presence of noise. When the number of componentsis known in advance
% Unlike DFT, it is able to estimate frequencies with accuracy higher than one sample, because its estimation function can be evaluated for any frequency, not just those of DFT bins. 
% This is a form of superresolution

% Its chief disadvantage is that it requires the number of components to be known in advance,
% so the original method cannot be used in more general cases. Methods exist for estimating 
% the number of source components purely from statistical properties of the autocorrelation matrix.
% See, e.g. [5] In addition, MUSIC assumes coexistent sources to be uncorrelated, which limits 
% its practical applications.


%%%%%%%%% -- Functions -- %%%%%%%%% 
function r_k = acf(x, biased)
% acf computes the estimate of a given input sequence. Gives user 
% choice between biased and unbiased estimates
    r_k = []; N = length(x);
    for k=0:N-1
        if strcmp(biased, 'biased')
            a = (1/N); % factor for biased estimate
        else
            a = (1/(N-k)); % factor for unbiased estimate
        end
        r_k(k+1) = a * dot(x(k+1:N), conj(x(1:N-k)));
    end
end
function [P1, f] = correlogram(x, biased, fs) 
    P1 = abs(fft(acf(x, biased))); L = length(x);
    P1 = P1(1:floor(L/2)+1);
    f = (fs/2)*[0:length(P1)- 1]/length(P1);
end