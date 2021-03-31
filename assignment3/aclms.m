function [q, error, y_hat] = aclms(y, x, M, lr, leak)
%ANC_LMS Implements ANC with LMS
%   Outputs:
%   q - dynamics of the weights learnt from the data
%   Inputs:
%   x - input signal
%   y - ouput signal
%   lr - step size/ learning rate for standard, initial for adaptive algs.

    q = zeros(2*(M+1), length(x),'like',1i);  
    error = zeros(size(x),'like',1i);
    y_hat = zeros(size(x),'like',1i);
    x_pad = [zeros(M,1,'like',1i); x]; % zero-padding
    for n = 1:length(x)
        x_a = [x_pad(n+M:-1:n); conj(x_pad(n+M:-1:n))];
        y_hat(n) = q(:,n)'*x_a;
        error(n) = y(n) - y_hat(n);
        if n < length(x)
            q(:, n+1) = (1-lr*leak)*q(:, n) + lr*conj(error(n))*x_a;
        end
    end
end

