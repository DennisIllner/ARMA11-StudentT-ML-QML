function [likelihood_contributions] = FUN_Likelihood_Students_t(parameterset,y)

% The Fun_Likelihood_Students_t function takes the parameterset (c,phi,theta,nu) of an
% ARMA(1,1) with student's t innovations  and the respective
% time series y as input values. 
% The function first backs out the series of innovations using the time
% series data y and assuming the first innovation e_1 to be 0. The first
% value then is dropped.
% Secondly the individual log likelihood contributions are calculated using the
% function of the student's t distribution. The function returns the
% individual log likelihood contributions.


%Extract Parameters from Parameterset
c     = parameterset(1);
phi   = parameterset(2);
theta = parameterset(3);
nu    = parameterset(4);

%Initialize vector for time series of innovations
e = zeros(length(y),1);

%Initialize fist innovation to be 0
%e(1) = 0 ;

%Loop to back out each innovation
for i=2:length(e);
    e(i)= y(i)-c-phi*y(i-1)-theta*e(i-1);
end

%Drop e_1 from vector of innovations
e =  e(2:end);


%Now Calculating the log-likelihood contribution:
likelihood_contributions = log(gamma((nu+1)/2))  - log(sqrt(nu*pi)*gamma(nu/2))  - ((nu+1)/2) * log(1 + ((e.^2)/nu));


end