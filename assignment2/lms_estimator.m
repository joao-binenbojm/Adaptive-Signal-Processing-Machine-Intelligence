function [params, error] = lms_estimator(data, order, lr)
%LMS_ESTIMATOR Computes the lms estimate using gradient descent
%   params - dynamics of the weights learnt from the data
%   error - dynamics of the prediction error
%   data - input containing samples of the process for estimation
%   order - order parameter that determines number of parameters to be
%   estimated
%   lr - step size/ learning rate 

params = zeros(order, length(data)); % n_parameters = order
error = zeros(size(data)) ;

for i = order+1:length(data)
    current_error = data(i) - dot(flip(data(i-order:i-1)), params(:, i-1));
    error(i) = current_error;
    params(:, i) = params(:, i-1) + lr*(current_error)*flip(data(i-order:i-1));
end

