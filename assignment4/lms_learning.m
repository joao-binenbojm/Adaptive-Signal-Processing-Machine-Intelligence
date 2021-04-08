function [params, err, y_hat] = lms_learning(y, x, order, lr, gamma, opt, lr_a, pt)
%LMS Computes the lms estimate using gradient descent
%   params - dynamics of the weights learnt from the data
%   err - dynamics of the prediction error
%   y_hat - predicted output by filter
%   y - desired signal
%   x - input signal used to generate output
%   order - order parameter that determines number of parameters to be
%   estimated
%   lr - step size/ learning rate 
%   gamma - leakage coefficient
%   opt - determines which lms algorithm to implement
%   lr_a - learning rate for scaling coefficient
%   pt - boolean determines whether to pretrain

    N = length(y);
    err = zeros(N,1);
    y_hat = zeros(N,1);

    % Pre-training weights
    if pt
        [params, a] = pretrain(y,x,100,20,order,lr,gamma,lr_a,opt);
    else
         params = zeros(order+1,N); % params + bias
         a = 30;
    end
    % Training on rest of data
    for n = 1:N
        if strcmp(opt,'standard')
            y_hat(n) = params(:,n).'*[1;x(n+order-1:-1:n)];
            act_factor = 1;
        elseif strcmp(opt,'nonlinear')
            nonlinear = tanh(params(:,n).'*[1;x(n+order-1:-1:n)]);
            y_hat(n) = a*nonlinear;
            act_factor = a*(1-(nonlinear.^2));
        else
            error("Please enter a valid algorithm option");
        end
        err(n) = y(n) - y_hat(n);
        if n < N
            params(:, n+1) = (1-lr*gamma)*params(:, n)...
                + lr*err(n)*act_factor*[1;x(n+order-1:-1:n)];
            if strcmp(opt,'nonlinear') % update scaling coefficient
                a = a + lr_a*err(n)*nonlinear;
            end
        end
    end
end

function [params, a] = pretrain(y,x,n_epochs,K,order,lr,gamma,lr_a,opt)
% this function is used to pretrain weights, bias and scaling coeff
    N = length(y);
    params = zeros(order+1,N); % params + bias
    weights = zeros(order+1,1);
    err = zeros(N,1);
    y_hat = zeros(N,1);
    a = 40; % initial scaling coefficient
    for e = 1:n_epochs
        for k = 1:K
            if strcmp(opt,'standard')
                y_hat(k) = weights.'*[1;x(k+order-1:-1:k)];
                act_factor = 1;
            elseif strcmp(opt,'nonlinear')
                nonlinear = tanh(weights.'*[1;x(k+order-1:-1:k)]);
                y_hat(k) = a*nonlinear;
                act_factor = a*(1-(nonlinear.^2));
            else
                error("Please enter a valid algorithm option");
            end
            err(k) = y(k) - y_hat(k);
            if k < K
                weights = (1-lr*gamma)*weights...
                    + lr*err(k)*act_factor*[1;x(k+order-1:-1:k)];
                if strcmp(opt,'nonlinear') % update scaling coefficient
                    a = a + lr_a*err(k)*nonlinear;
                end
            end
        end
    end
    params(:,1:K) = repmat(weights, 1, K); % setting best lerned params 
end
