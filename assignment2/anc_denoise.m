function xhat = anc_denoise(signal, sec_noise, w)
%ANC_DENOISE: Implements ANC with LMS
%   Outputs:
%   xhat - denoised signal
%   Inputs:
%   signal - noise-corrupted signal
%   sec_noise - secondary noise input (reference noise input)
%   w - learned weights between secondary and primary noise

    eta = zeros(size(signal));
    xhat = zeros(size(signal));
    M = length(w);
    u = delayseq(repmat(sec_noise, 1, M), [0:M-1])';
    
    for n = 1:length(signal)
        eta(n) = dot(w, u(:, n));
        xhat(n) = signal(n)  - eta(n);
    end
end

