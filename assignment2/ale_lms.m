function [w, xhat, error] = ale_lms(signal, lr, delta, M)
%ALE_LMS Implements ALE with LMS
%   Outputs:
%   w - dynamics of the weights learnt from the data
%   error - dynamics of the prediction error
%   Inputs:
%   signal - noise-corrupted signal
%   lr - step size/ learning rate for standard, initial for adaptive algs.
%   delta - delay between current signal and predictor input
%   M - ALE filter length

    w = zeros(M, length(signal)+1); 
    error = zeros(size(signal));
    xhat = signal(size(signal));
    
    for n = delta+M:length(signal)
        % Define current u-vector
        u = flip(signal(n-delta-M+1:n-delta));
        xhat(n) = dot(w(:, n), u);
        error(n) = signal(n)  - xhat(n);
        % Update weights
        w(:, n+1) = w(:, n) + lr*error(n)*u;
    end
end

