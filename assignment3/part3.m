%% --- 3. Adaptive signal processing --- %%

%% 3.3. A Real Time Spectrum Analyser Using Least Mean Square %%
%(a)
clc; clear variables; close all;

n_samples = 1500;
f = [100*ones(500,1); ...
    100*ones(500,1) + ([501:1000]' - 500)/2; ...
    100*ones(500,1) + (([1001:1500]' -  1000)/25).^2];
phi = cumsum(f);
fs = 1000;
variance = 0.05;
eta = wgn(n_samples,1,pow2db(variance),'complex');
y = exp(1i*2*pi*phi/fs) + eta;
% Creating fourier signal
time = 1:n_samples;
f_ax = 0:fs-1;
x = (1/fs)*exp(1i*2*pi*f_ax'*time./fs);

% Testing out and plotting for multiple different leaks
leaks = [0, 0.25, 0.5];
figure(1);
for idx = 1:3
    [w,~] = dft_clms(y.', x, leaks(idx));
    % Plotting results
    subplot(3,1,idx); hold on; set(gca,'fontsize', 16);
    mesh(time, f_ax(1:floor(fs/2)), abs(w(1:floor(fs/2),:)).^2);
    view(2);
    xlabel('Time (samples)');
    ylabel('Frequency (Hz)');
    ylim([0,500]);
    title(strcat('CLMS Fourier Estimate ($\gamma=',num2str(leaks(idx)),'$)'),'Interpreter','Latex');
    hold off;
end

%% (d)
clc; clear variables; close all;
load('EEG_Data_Assignment2.mat');

a = 3000;
t_range = a:a+1199;
POz = detrend(POz(t_range));
n_samples = length(t_range);
f_ax = 0:fs-1;
x = (1/fs)*exp(1i*2*pi*f_ax'*t_range./fs);
[w,~] = dft_clms(POz, x, 0);
% Plotting results
figure(1); hold on; set(gca,'fontsize', 16);
mesh(t_range, f_ax(1:floor(fs/2)), abs(w(1:floor(fs/2),:)).^2);
view(2); ylim([0,60]);
xlabel('Time (samples)');
ylabel('Frequency (Hz)');
title('CLMS Fourier Estimate (POz)');

  