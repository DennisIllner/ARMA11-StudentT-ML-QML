function [likelihood_contributions] = FUN_QML_Gaussian(parameterset,y)

% The Fun_Likelihood_Gaussian function takes a parameter set (c, phi, theta, sigma^2) 
% and the corresponding time series y as input. 
%
% First, the function extracts the series of innovations using the time 
% series data y, assuming the first innovation e_1 to be 0. The first value 
% is then dropped from the innovation series.
%
% Second, the function calculates the individual log-likelihood contributions 
% using the probability density function of the normal distribution. The function 
% returns the individual log-likelihood contributions.
%
% Note that QML (Quasi-Maximum Likelihood) is used when the true underlying 
% distribution of the stochastic process is unknown. In this context, we assume 
% that the normal distribution could be a plausible candidate, though we are not 
% certain that it is the correct one.

%Extract Parameters from Parameterset
c     = parameterset(1);
phi   = parameterset(2);
theta = parameterset(3);
sigma = parameterset(4);

%Initialize vector for time series of innovations
e = zeros(length(y),1);

%Initialize fist innovation to be 0
e(1) = 0 ;

%Loop to back out each innovation
for i=2:length(e);
    e(i)= y(i)-c-phi*y(i-1)-theta*e(i-1);
end

%Drop e_1 from vector of innovations
e =  e(2:end);


%Now Calculating the log-likelihood contribution:
likelihood_contributions = log(1/sqrt(2*pi*sigma))-((e.^2)/(2*sigma));


end