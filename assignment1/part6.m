%% --- 1. CLASSICAL AND MODERN SPECTRUM ESTIMATION --- %%

%% 1.5 Real World Signals: Respiratory Sinus Arrhythmia from RR-Intervals %%

%% Part (a)
clc; clear all; close all;
load PCAPCR;
[Ux, Sx, Vx] = svd(X); [Un, Sn, Vn] = svd(Xnoise);

% Plotting singular values of X and Xnoise
figure(1); subplot(1, 2, 1); hold on;
stem(diag(Sx)); stem(diag(Sn)); set(gca,'fontsize', 14);
xlabel('Subspace Dimension Index'); ylabel('Singular Value');
title('Singular Values of Dataset'); legend('X', 'Xnoise'); hold off;

% Square error between SVs
subplot(1, 2, 2);
stem((diag(Sx) - diag(Sn)).^2); set(gca,'fontsize', 14);
xlabel('Subspace Dimension Index'); ylabel('Singular Value Square Error');
title('Singular Value Square Errors');

find_components(diag(Sn), 0.9)
find_components(diag(Sx), 0.9)
find_components(diag(Sn), 0.95)
find_components(diag(Sx), 0.95)

%% Part (b)
clc; clear all; close all;
load PCAPCR;

[U, S, V] = svds(Xnoise, 3); % keeps only principal components (3)
X_noise_approx = U*S*V';

% Plotting MSE per column or noise and noise_approx
figure(1); hold on;
stem(mean(((X - Xnoise).^2), 1)); stem(mean(((X - X_noise_approx).^2), 1)); 
set(gca,'fontsize', 14); xlabel('Columns'); ylabel('MSE')
title('MSE Between Original and Noisy Data');
legend('X, X$_{noise}$', 'X, $\tilde{X}_{noise}$', 'Interpreter','latex'); hold off;
immse(X, Xnoise)
immse(X, X_noise_approx)

    %% Part (c)
clc; clear all; close all;
load PCAPCR;

% Computing weights
B_ols = inv(Xnoise' * Xnoise)*X'*Y;
B_pcr = pcr(Y, Xnoise, 3);

% Testing results
Y_ols = Xnoise * B_ols;
Y_pcr = Xnoise * B_pcr;

Ytest_ols = Xtest * B_ols;
Ytest_pcr = Xtest * B_pcr;

% Plotting singular values of X and Xnoise
figure(1); subplot(1, 2, 1); hold on;
stem(sum((Y-Y_ols).^2, 1)); stem(sum((Y-Y_pcr).^2, 1)); set(gca,'fontsize', 14);
xlabel('Outputs'); ylabel('Squared Errors');
title('Estimation Error'); legend('OLS', 'PCR'); hold off;

% Square error between SVs
subplot(1, 2, 2); hold on;
stem(sum((Ytest-Ytest_ols).^2, 1)); stem(sum((Ytest-Ytest_pcr).^2, 1)); set(gca,'fontsize', 14);
xlabel('Outputs'); ylabel('Squared Errors');
title('Test Error'); legend('OLS', 'PCR'); hold off;

fprintf('OLS Metrics (Estimation): \n RMSE: %.3f \n MAE: %.3f \n\n', myrmse(Y, Y_ols), mymae(Y, Y_ols))
fprintf('PCR Metrics (Estimation): \n RMSE: %.3f \n MAE: %.3f \n\n', myrmse(Y, Y_pcr), mymae(Y, Y_pcr))

fprintf('OLS Metrics (Test): \n RMSE: %.3f \n MAE: %.3f \n\n', myrmse(Ytest, Ytest_ols), mymae(Ytest, Ytest_ols))
fprintf('PCR Metrics (Test): \n RMSE: %.3f \n MAE: %.3f \n\n', myrmse(Ytest, Ytest_pcr), mymae(Ytest, Ytest_pcr))


% Part (d)
ols_error = [];
pcr_error = [];

for i = 1:100
    [Y_ols, Y_ols_hat] = regval(B_ols);
    [Y_pcr, Y_pcr_hat] = regval(B_pcr);
    ols_error = [ols_error; sum((Y_ols - Y_ols_hat).^2, 1)];
    pcr_error = [pcr_error; sum((Y_pcr - Y_pcr_hat).^2, 1)];
end

ols_error = mean(ols_error, 1);
pcr_error = mean(pcr_error, 1);

figure(2); hold on;
stem(ols_error); stem(pcr_error); 
set(gca,'fontsize', 14); xlabel('Columns'); ylabel('Averaged Square Error')
title('Averaged Square Error Between Estimate & Original Data');
legend('OLS', 'PCR'); hold off;



function B = pcr(Y, X, k)
    [U, S, V] = svds(X, 3); % only keep r largest components
    for i = 1:length(diag(S)) % Computing pseudoinverse
        S(i, i) = 1/S(i, i);
    end
    B = (V * S' * U') * Y;
end
function rmse = myrmse(X1, X2)
    rmse = mean((X1 - X2).^2, 'all').^0.5;
end
function mae = mymae(X1, X2)
    mae = mean(abs(X1 - X2), 'all');
end
function n_components = find_components(SVs, proportion) % Find number of
    eig_var = sum(SVs); 
    i = 1;
    while sum(SVs(1:i))/eig_var < proportion
        i = i + 1;
        n_components = i;
    end
end