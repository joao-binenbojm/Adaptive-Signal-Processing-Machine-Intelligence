function [ar_params, ma_params, error] = lms_arma(output, input, p, q, lr, gass, rho, alpha, leak)
%LMS_ESTIMATOR Computes the lms estimate using gradient descent
%   Outputs:
%   ar_params - dynamics of the AR weights learnt from the data
%   ar_params - dynamics of the MA weights learnt from the data
%   error - dynamics of the prediction error
%   Inputs:
%   output - contains samples of output process
%   input - contains samples of input process
%   p - AR order parameter
%   q - MA order parameter
%   lr - step size/ learning rate for standard, initial for adaptive algs.
%   gass - determines what type of adaptive LMS weight update is carried
%       'standard' computes the standard non-adaptive lms algorithm
%       'ben' computes the Beneviste update
%       'af' computes the Ang & Farhang update
%       'mx' computes the Matthews & Xie update
%   rho - optional parameter if non-standard updates are chosen
%   alpha - optional parameter if 'a&f' update is chosen
%   leak - leakage coefficient


    params = zeros(p+q, length(output)); % n_parameters = order
    phi = zeros(p+q, length(output));
    lrs = lr*ones(size(output));
    error = ones(size(output));
    
    for i = max([p,q])+1:length(output)-1
        % Create augmented data vector with both y and x values
        aug_dat = [flip(output(i-p:i-1)), flip(input(i-q:i-1))]';
        error(i) = output(i) - dot(aug_dat, params(:, i));
        % SimultaneousAR and MA updates
        params(:, i+1) = (1-leak*lr)*params(:, i) + lrs(i)*(error(i))*aug_dat;
        % Step-size update
        if strcmp(gass, 'af')
            if i > max([p,q])+1
                phi(:, i) = alpha*phi(:, i-1) + error(i-1)*prev_aug_dat;
            end
            lrs(i+1) = lrs(i) + rho*error(i)*dot(aug_dat, phi(:, i));
        elseif strcmp(gass, 'ben')
            if i > max([p,q])+1
                phi(:, i) = (eye(length(prev_aug_dat))-lrs(i-1)*(prev_aug_dat(:)*prev_aug_dat(:).'))...
                    *phi(:, i-1) + error(i-1)*prev_aug_dat;
            end
            lrs(i+1) = lrs(i) + rho*error(i)*dot(aug_dat, phi(:, i));
        elseif strcmp(gass, 'mx')
            if i > max([p,q])+1
                phi(:, i) = error(i-1)*prev_aug_dat;
            end
            lrs(i+1) = lrs(i) + rho*error(i)*dot(aug_dat, phi(:, i));
        elseif strcmp(gass, 'standard')
        else
            error('Either wrong number of inputs or wrong step update method selected');
        end
        prev_aug_dat = aug_dat;
    end
    ar_params = params(1:p, :); ma_params = params(p+1:q, :);
end

