%% --- 3. Adaptive signal processing --- %%

%% 3.2. Adaptive AR Model Based Time-Frequency Estimation %%

%% (a) 
clc; clear variables; close all;

n_samples = 1500;
f = [100*ones(500,1); ...
    100*ones(500,1) + ([501:1000]' - 500)/2; ...
    100*ones(500,1) + (([1001:1500]' -  1000)/25).^2];
phi = cumsum(f);
fs = 1500;
variance = 0.05;
eta = wgn(n_samples,1,pow2db(variance),'complex');
y = exp(1i*2*pi*phi/fs) + eta;
% Plotting frequency
figure(1); subplot(1,2,1); hold on; set(gca,'fontsize', 16);
plot([1:n_samples], f);
xlabel('Time (samples)'); ylabel('Frequency (Hz)');
title('Frequency modulated (FM) signal'); hold off;
% Plotting PSD of signal from estimates
order = [1, 5, 10, 20];
subplot(1,2,2); hold on; set(gca,'fontsize', 16);
H = zeros(4,1500); w = zeros(4,1500);
for i = 1:4
    % Finding AR coefficients
    a = aryule(y,order(i));
    [H(i,:),w(i,:)] = freqz(1,a,n_samples);
    P = abs(H(i,:)).^2;
    plot(w(i,:)/pi, pow2db(P));
end
title('AR Spectrum estimates')
xlabel('Normalized Frequency (\times \pi rad/sample)'); 
ylabel('Power/frequency (dB/(rad/sample))');
legend(arrayfun(@(a)strcat('p=',num2str(a)),order,'uni',0)); 
hold off;
% Print the mean frequency of signal and mean freq of peak
fprintf(mean(f));
[~,I] = max(H(1,:));
fprintf((w(1,I)/pi)*(fs/2));

%% (b)
clc; clear variables; close all;

n_samples = 1500;
f = [100*ones(500,1); ...
    100*ones(500,1) + ([501:1000]' - 500)/2; ...
    100*ones(500,1) + (([1001:1500]' -  1000)/25).^2];
phi = cumsum(f);
fs = 1500;
variance = 0.05;
eta = wgn(n_samples,1,pow2db(variance),'complex');
y = exp(1i*2*pi*phi/fs) + eta;

% CLMS coefficients
x = delayseq(y,1);
lrs = [1e-1, 1e-2,1e-3];
figure(1);
for idx = 1:3
    [a,~,~] = clms(y, x, 0, lrs(idx), 0);
    H = zeros(n_samples,n_samples);
    for n = 1:n_samples
        % Run complex-valued LMS algorithm to estimate AR coefficient aˆ1(n)
        [h, w] = freqz(1 , [1; -conj(a(n))], n_samples); % Compute power spectrum
        H(:, n) = abs(h).^2; % Store it in a matrix 
    end
    % Remove outliers in the matrix H
    medianH = 50*median(median(H));
    H(H > medianH) = medianH;
    % Plotting
    subplot(1,3,idx); hold on; set(gca,'fontsize', 16);
    mesh([1:n_samples], w/pi, H);
    view(2);
    xlabel('Time (samples)');
    ylabel('Frequency (\times \pi rad/sample)');
    title(strcat('CLMS Spectral Estimate ($\mu=',num2str(lrs(idx)),'$)'), 'Interpreter', 'Latex');
end

