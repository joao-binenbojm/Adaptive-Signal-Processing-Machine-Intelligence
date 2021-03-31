function [eta,rho] = circularity(data)
%CIRCULARITY(data) Computes circularity quotient and coefficient
%   Takes in complex data and returns circularity coefficient and
%   circularity quotient
    p = mean(data.*data); c = mean(abs(data).^2);
    rho = p/c;
    eta = abs(rho);
end

