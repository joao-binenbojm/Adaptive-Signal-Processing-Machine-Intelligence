function [w, xhat] = anc_lms(signal, sec_noise, lr, M)
%ANC_LMS Implements ANC with LMS
%   Outputs:
%   w - dynamics of the weights learnt from the data
%   xhat - denoised signal dynamics
%   Inputs:
%   signal - noise-corrupted signal
%   sec_noise - secondary noise input (reference noise input)
%   lr - step size/ learning rate for standard, initial for adaptive algs.
%   M - ANC secondary noise input length

    w = zeros(M, length(signal)); 
    eta = zeros(size(signal));
    xhat = zeros(size(signal));
    u = delayseq(repmat(sec_noise, 1, M), [0:M-1])';
 
    for n = 1:length(signal)
        eta(n) = dot(w(:, n), u(:, n));
        xhat(n) = signal(n)  - eta(n);
        % Update weights
        w(:, n+1) = w(:, n) + lr*xhat(n)*u(:, n);
    end
end

