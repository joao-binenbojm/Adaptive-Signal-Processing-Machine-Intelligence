function [ar_params, ma_params, error] = gngd(output, input, p, q, lr, rho, leak)
%GNGD Computes the GNGD estimate using gradient descent
%   Outputs:
%   ar_params - dynamics of the AR weights learnt from the data
%   ar_params - dynamics of the MA weights learnt from the data
%   error - dynamics of the prediction error
%   Inputs:
%   output - contains samples of output process
%   input - contains samples of input process
%   p - AR order parameter
%   q - MA order parameter
%   lr - step size/ learning rate 
%   rho - step size for regularization updates
%   leak - leakage coefficient

    params = zeros(p+q, length(output)); % n_parameters = order
    error = ones(size(output));
    reg = ones(size(output))/lr;
    
    for i = max([p,q])+1:length(output)-1
        % Create augmented data vector with both y and x values
        aug_dat = [flip(output(i-p:i-1)), flip(input(i-q:i-1))]';
        error(i) = output(i) - dot(aug_dat, params(:, i));
        % Simultaneous AR and MA updates
        lr_now = lr/(reg(i) + dot(aug_dat, aug_dat));
        params(:, i+1) = (1-leak*lr)*params(:, i) + lr_now*(error(i))*aug_dat;
        % Update regularization parameter
        if i >  max([p,q])+1
            num = rho*lr*error(i)*error(i-1)*dot(old_dat, aug_dat);
            den = ( reg(i-1) + dot(old_dat, old_dat) ).^2;
            reg(i+1) = reg(i) - num/den;
        end
        old_dat = aug_dat;
        
    end
    ar_params = params(1:p, :); ma_params = params(p+1:q, :);
end

