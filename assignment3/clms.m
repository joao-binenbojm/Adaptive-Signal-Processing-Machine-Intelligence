function [h, error, y_hat] = clms(y, x, M, lr, leak)
%ANC_LMS Implements ANC with LMS
%   Outputs:
%   h - dynamics of the weights learnt from the data
%   y_hat - denoised signal dynamics
%   Inputs:
%   x - noise-corrupted signal
%   lr - step size/ learning rate for standard, initial for adaptive algs.

    h = zeros(M+1, length(x),'like',1i); 
    error = zeros(size(x),'like',1i);
    y_hat = zeros(size(x),'like',1i);
    x_pad = [zeros(M,1); x]; % zero-padding
    for n = 1:length(x)
        y_hat(n) = h(:,n)'*x_pad(n+M:-1:n); 
        error(n) = y(n) - y_hat(n);
        if n < length(x)
            h(:, n+1) = (1-lr*leak)*h(:, n) + lr*conj(error(n))*x_pad(n+M:-1:n);
        end
    end
end

