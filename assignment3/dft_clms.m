function [h, error] = dft_clms(y, x, leak)
%DFT_CLMS(y,x,M,lr,leak) daptive spectral estimte using CLMS
%   Outputs:
%   h - dynamics of the weights learnt from the data
%   y_hat - denoised signal dynamics
%   Inputs:
%   x - Fourier input
%   y - ouput signal
%   lr - step size/ learning rate for standard, initial for adaptive algs.
%   leak - degree of leakage during learning

    [nft, n_samples] = size(x);
    h = zeros(nft, n_samples,'like',1i); 
    error = zeros(n_samples,1,'like',1i);
    lr = 1;
    for n = 1:n_samples
        error(n) = y(n) - h(:,n)' * x(:,n);
        if n < length(x)
            h(:, n+1) = (1-lr*leak)*h(:, n) + lr*conj(error(n))*x(:,n);
        end
    end
end

