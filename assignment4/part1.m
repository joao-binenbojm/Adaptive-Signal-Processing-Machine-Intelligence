%% --- 4. From LMS to Deep Learning --- %%
%% (1) 
clc; clear variables; close all;

% Add data folder to path
addpath('../data');
load('time-series.mat');

% Learning parameters
y = detrend(y); % remove mean
order = 4;
x = [zeros(order,1); y]; % input to AR(4) model
[params,error,y_hat] = lms_learning(y, x, order, 1e-5, 0, 'standard',0,false);

figure(1); 
subplot(1,2,1); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y);
plot(1:length(y_hat), y_hat);
xlabel('Samples (n)');
ylabel('Y'); title('AR(4) Prediction');
legend('True Signal', 'AR(4)'); hold off;
subplot(1,2,2); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y);
plot(1:length(y_hat), y_hat);
xlim([800,1000]);
xlabel('Samples (n)');
ylabel('Y'); title('AR(4) Prediction (Segment)');
legend('True Signal', 'AR(4)'); hold off;

% Computing MSE and prediction gain
MSE = (error'*error)/1000;
Rp = pow2db(var(y_hat)/var(error));
fprintf(sprintf('MSE: %.4f, Rp: %.4f',MSE,Rp));

%% (2/3) 
clc; clear variables; close all;

% Add data folder to path
addpath('../data');
load('time-series.mat');

% Learning parameters
y = detrend(y); % remove mean
order = 4;
x = [zeros(order,1); y]; % input to AR(4) model
[params,error,y_hat] = lms_learning(y, x, order, 2e-7, 0, 'nonlinear',0.1,false);

figure(1); 
subplot(1,2,1); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y);
plot(1:length(y_hat), y_hat);
xlabel('Samples (n)');
ylabel('Y'); title('Scaled LMS-Tanh Prediction');
legend('True Signal', 'LMS + tanh'); hold off;
subplot(1,2,2); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y);
plot(1:length(y_hat), y_hat);
xlim([800,1000]);
xlabel('Samples (n)');
ylabel('Y'); title(' Scaled LMS-Tanh Prediction (Segment)');
legend('True Signal', 'LMS + tanh'); hold off;

% Computing MSE and prediction gain
MSE = (error'*error)/1000;
Rp = pow2db(var(y_hat)/var(error));
fprintf(sprintf('MSE: %.4f, Rp: %.4f',MSE,Rp));

%% (4) adding bias
clc; clear variables; close all;

% Add data folder to path
addpath('../data');
load('time-series.mat');

% Learning parameters
order = 4;
x = [zeros(order,1); y]; % input to AR(4) model
[params,error,y_hat] = lms_learning(y, x, order, 2e-7, 0, 'nonlinear',1e-1,false);

figure(1); 
subplot(1,2,1); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y);
plot(1:length(y_hat), y_hat);
xlabel('Samples (n)');
ylabel('Y'); title('Prediction with bias');
legend('True Signal', 'LMS + tanh'); hold off;
subplot(1,2,2); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y);
plot(1:length(y_hat), y_hat);
xlim([800,1000]);
xlabel('Samples (n)');
ylabel('Y'); title('Prediction with bias (Segment)');
legend('True Signal', 'LMS + tanh'); hold off;

% Computing MSE and prediction gain
MSE = (error'*error)/1000;
Rp = pow2db(var(y_hat)/var(error));
fprintf(sprintf('MSE: %.4f, Rp: %.4f',MSE,Rp));

%% (5) Introducted pretraining
clc; clear variables; close all;

% Add data folder to path
addpath('../data');
load('time-series.mat');

% Learning parameters
order = 4;
x = [zeros(order,1); y]; % input to AR(4) model
[params,error,y_hat] = lms_learning(y, x, order, 2e-7, 0, 'nonlinear',1e-1,true);

figure(1); 
subplot(1,2,1); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y);
plot(1:length(y_hat), y_hat);
xlabel('Samples (n)');
ylabel('Y'); title('Prediction with pretraining');
legend('True Signal', 'LMS + tanh'); hold off;
subplot(1,2,2); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y);
plot(1:length(y_hat), y_hat);
xlim([800,1000]);
xlabel('Samples (n)');
ylabel('Y'); title('Prediction with pretraining (Segment)');
legend('True Signal', 'LMS + tanh'); hold off;

% Computing MSE and prediction gain
MSE = (error'*error)/1000;
Rp = pow2db(var(y_hat)/var(error));
fprintf(sprintf('MSE: %.4f, Rp: %.4f',MSE,Rp));

%% (6) see report

%% (7) see report

%% (8) see report