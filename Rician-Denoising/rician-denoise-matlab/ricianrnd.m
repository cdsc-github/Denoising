function nu = ricianrnd(nu,sigma)
%RICIANRND  Random arrays from Rician distribution.
%   r = RICIANRND(nu,sigma) returns random numbers from the Rician
%   distribution with parameters nu and sigma.  If nu is an array, then r
%   is an array of the same shape.

% Pascal Getreuer 2009


r = nu(:)'/sqrt(2);
r = (r([1,1],:) + randn(2,length(r))*sigma).^2;
nu(:) = sqrt(sum(r));
