%% --- 3. Adaptive signal processing --- %%

%% 3.1. Complex LMS and Widely Linear Modelling %%

%% (a) 
clc; clear variables; close all;
b1 = 1.5 + 1i; b2 = 2.5 - 0.5i;
realizations = 100;
input = wgn(1000,realizations,0,'complex');
data = zeros(1000,realizations,'like',1i);
error_store = zeros(2, 1000,realizations);
f = waitbar(0,'Processing...');
for idx = 1:realizations
    for jdx = 1:1000
        if jdx == 1
            data(jdx,idx) = 0i;
        else
            data(jdx,idx) = input(jdx,idx) + b1*input(jdx-1,idx) + b2*conj(input(jdx-1,idx));
        end
    end
    [~,error_store(1, :, idx),~] = aclms(data(:,idx),input(:,idx),1, 0.1, 0);
    [~,error_store(2, :, idx),~] = clms(data(:,idx),input(:,idx),1, 0.1, 0);
    waitbar(idx/realizations,f,'Processing...');
end
close(f);
mean_error = mean(abs(error_store).^2,3);
% Plots
figure(1); subplot(1,3,1); hold on; set(gca,'fontsize', 18);
plot(input, 'bo'); title('WGN'); 
xlabel('Re(z)'); ylabel('Im(z)'); hold off;
subplot(1,3,2); hold on; set(gca,'fontsize', 18);
plot(data, 'ro'); title('WLMA(1)'); 
xlabel('Re(z)'); ylabel('Im(z)'); hold off;
subplot(1,3,3); hold on; set(gca,'fontsize', 18);
plot([1:1000], pow2db(mean_error(1,:)));
plot([1:1000], pow2db(mean_error(2,:)));
ylabel('Squared Error (dB)'); xlabel('Time (samples)');
title('Learning Curve'); legend('ACMLS', 'CMLS');
hold off;

%% (b)
clc; clear variables; close all;

% Combine two components as complex numbers and plot scatter plots
complex_v = zeros(3,5000,'like',1i);
names = {'high-wind','medium-wind','low-wind'};
markers = ["or", "ob", "ok"];
figure(1);
for i = 1:3
    % Loading data
    data = load(names{i});
    complex_v(i,:) = complex(data.v_east,data.v_north);
    % Parsing data strings for title
    str = split(names{i},'-');
    tit = strcat(upper(str{1}(1)), str{1}(2:end), ' Wind');
    % Circularity plots
    subplot(1,3,i); hold on; set(gca,'fontsize', 18);
    plot(complex_v(i,:), markers(i));
    xlabel('Re(z)'); ylabel('Im(z)');
    title(tit);
    hold off;
end

% Grid search to find optimal filter length for both algorithms on each
% data set
% Computing values and creating plots
n_orders = 20;
MSPE = zeros(2,3,n_orders); 
lrs = [0.001, 0.01,0.1]; % for a given wind type
counter = 0; % subplot counter
figure(2);
for j = 1:3
    counter = counter + 1 ;
    subplot(1, 3,counter); hold on; set(gca,'fontsize', 18);
    for i = 1:2
        input = delayseq(complex_v(j,:).',1);
        for l = 1:n_orders
            if i == 1
                [~,err,~] = ...
                    clms(complex_v(j,:).',input,l-1, lrs(j), 0);
            else
                [~,err,~] = ...
                    aclms(complex_v(j,:).',input,l-1, lrs(j), 0);
            end
            sq_err = abs(err).^2;
            MSPE(i,j,l) = mean(sq_err);
        end
        plot(1:n_orders, 10*log10(squeeze(MSPE(i,j,:))));
        % Parsing data strings for title
    end
    str = split(names{j},'-');
    tit = strcat(upper(str{1}(1)), str{1}(2:end), ' Wind');
    legend('CLMS','ALMS');
    xlabel('Filter Order'); ylabel('MSPE (dB)');
    title(tit);
    hold off;
end

%% (c)
clc; clear variables; close all;

% Assume zero distortions, but mismatch in magnitude

f0 = 50; % nominal frequency
fs = 5000; % sampling frequency
V = 10; % balanced magnitude
phi = pi/3; % phase shift
% Generate balanced complex voltages
n = [1:10000]; % one second of data in this made up system
v = sqrt(3/2)*V*exp(i*2*pi*(f0/fs)*n + phi);
% Cicularity with changing Va and delta_a (20% jumps)
Va = [10:1:22]; delta_b = [0:pi/36:pi/3];
eta = zeros(2,length(Va));
for idx = 1:length(Va)
    % Generate unbalanced complex voltage 1
    A = (sqrt(6)/6)*(Va(idx) + 2*V);
    B = (sqrt(6)/6)*(Va(idx) + V*exp(-i*2*pi/3) + V*exp(i*2*pi/3));
    v_unbalanced1 = A*exp(i*2*pi*(f0/fs)*n + phi) + ...
    B*exp(-i*2*pi*(f0/fs)*n + phi);
    [eta(1,idx),~] = circularity(v_unbalanced1);
    % Generate unbalanced complex voltage 2
    A = (sqrt(6)/6)*V*(2 + exp(i*delta_b(idx)));
    B = (sqrt(6)/6)*V*(1 + exp(-i*(delta_b(idx) + 2*pi/3)) + exp(i*2*pi/3));
    v_unbalanced2 = A*exp(i*2*pi*(f0/fs)*n + phi) + ...
    B*exp(-i*2*pi*(f0/fs)*n + phi);
    [eta(2,idx),~] = circularity(v_unbalanced2);
end
% Plotting circularity diagrams
figure(1); subplot(1,2,1); 
hold on; set(gca,'fontsize', 16);
plot(v,'bo'); plot(v_unbalanced1,'ro'); plot(v_unbalanced2,'ko');
xlabel('Re(z)'); ylabel('Im(z)');
title('System Voltages');
legend('Balanced', 'Unbalanced (V)', 'Unbalanced ($\phi$)','Interpreter','Latex');
subplot(1,2,2); hold on; set(gca,'fontsize', 16);
plot([0:0.1:1.2], eta(1,:))
plot([0:0.1:1.2], eta(2,:))
xlabel('Deviation from balanced system (%)'); ylabel('Circularity Coefficient');
title('Divergence from balance and circularity');
legend('V','$\phi$','Interpreter','Latex');

%% (d)
clc; clear variables; close all;
% Proved like a B.O.S.S.

%% (e)
clc; clear variables; close all;
f0 = 50; % nominal frequency
fs = 500; % sampling frequency
V = 10; % balanced magnitude
phi = pi/3; % phase shift
% Generate balanced complex voltages
n = [1:200]; % one second of data in this made up system
v = sqrt(3/2)*V*exp(i*2*pi*(f0/fs)*n + phi);
% Unbalanced voltage params
Va = 20; delta_b = pi/3;
% Generate unbalanced complex voltage 1
A = (sqrt(6)/6)*(Va + 2*V);
B = (sqrt(6)/6)*(Va + V*exp(-i*2*pi/3) + V*exp(i*2*pi/3));
v_unbalanced1 = A*exp(i*2*pi*(f0/fs)*n + phi) + ...
B*exp(-i*2*pi*(f0/fs)*n + phi);
[eta(1),~] = circularity(v_unbalanced1);
% Generate unbalanced complex voltage 2
A = (sqrt(6)/6)*V*(2 + exp(i*delta_b));
B = (sqrt(6)/6)*V*(1 + exp(-i*(delta_b + 2*pi/3)) + exp(i*2*pi/3));
v_unbalanced2 = A*exp(i*2*pi*(f0/fs)*n + phi) + ...
B*exp(-i*2*pi*(f0/fs)*n + phi);
[eta(2),~] = circularity(v_unbalanced2);
% Estimation of different signals with CLMS and ACLMS
balanced_input = delayseq(v',1); % balanced system
[h_clms_balanced, err_clms_balanced, ~] = ...
    clms(v',balanced_input,0,0.0001,0);
[q_aclms_balanced, err_aclms_balanced, ~] = ...
    aclms(v',balanced_input,0,0.0001,0);
unbalanced_input1 = delayseq(v_unbalanced1', 1); % unbalanced system(V)
[h_clms_unbalanced1, err_clms_unbalanced1, ~] = ...
    clms(v_unbalanced1',unbalanced_input1,0,0.0001,0);
[q_aclms_unbalanced1, err_aclms_unbalanced1, ~] = ...
    aclms(v_unbalanced1',unbalanced_input1,0,0.00001,0);
unbalanced_input2 = delayseq(v_unbalanced2', 1); % unbalanced system(phase)
[h_clms_unbalanced2, err_clms_unbalanced2, ~] = ...
    clms(v_unbalanced2',unbalanced_input2,0,0.0001,0);
[q_aclms_unbalanced2, err_aclms_unbalanced2, ~] = ...
    aclms(v_unbalanced2',unbalanced_input2,0,0.0001,0);
% Frequency estimation
f_clms_balanced = FSL(h_clms_balanced,fs);
f_aclms_balanced = FWL(q_aclms_balanced,fs);
f_clms_unbalanced1 = FSL(h_clms_unbalanced1,fs);
f_aclms_unbalanced1 = FWL(q_aclms_unbalanced1,fs);
f_clms_unbalanced2 = FSL(h_clms_unbalanced2,fs);
f_aclms_unbalanced2 = FWL(q_aclms_unbalanced2,fs);
%Plotting frequency estimates
figure(1); subplot(1,3,1); hold on; set(gca,'fontsize', 16);
plot([1:length(f_clms_balanced)-2], f_clms_balanced(3:end),'LineWidth',3);
plot([1:length(f_aclms_balanced)-2],f_aclms_balanced(3:end),'LineWidth',3);
yticks([0,25,50,75,100]);
xlabel('Time (samples)'); ylabel('Frequency (Hz)');
title('Balaced System','Interpreter','Latex'); 
legend('CLMS','ACLMS'); ylim([0,100]);

subplot(1,3,2); hold on; set(gca,'fontsize', 16);
plot([1:length(f_clms_unbalanced1)-2], f_clms_unbalanced1(3:end),'LineWidth',2);
plot([1:length(f_aclms_unbalanced1)-2],f_aclms_unbalanced1(3:end),'LineWidth',2);
yticks([0,25,50,75,100]);
xlabel('Time (samples)'); ylabel('Frequency (Hz)');
title('Unbalanced System (V)','Interpreter','Latex'); 
legend('CLMS','ACLMS'); ylim([0,100]);

subplot(1,3,3); hold on; set(gca,'fontsize', 16);
plot([1:length(f_clms_unbalanced2)-2], f_clms_unbalanced2(3:end),'LineWidth',2);
plot([1:length(f_aclms_unbalanced2)-2],f_aclms_unbalanced2(3:end),'LineWidth',2);
yticks([0,25,50,75,100]);
xlabel('Time (samples)'); ylabel('Frequency (Hz)');
title('Unbalanced System ($\phi$)','Interpreter','Latex'); 
legend('CLMS','ACLMS'); ylim([0,100]);


% Frequency estimation of stricly linear autoregressive model
function f = FSL(h,fs)
    f = (fs/(2*pi))*atan(imag(h)./real(h));
end
% Frequency estimation of widely linear autoregressive model
function f = FWL(q,fs)
    h = q(1,:); 
    g = q(2,:);
    f = (fs/(2*pi))*...
        atan((sqrt(imag(h).^2 - abs(g).^2)./real(h)));
end