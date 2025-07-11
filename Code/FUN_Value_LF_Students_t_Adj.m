function [Value_LF] = FUN_Value_LF_Students_t_Adj(parameterset,y)
% The FUN_Value_LF_Students_t function calculates the value of the
% Likelihood function for the student's t distribution.
% The function takes the parameterset of an ARMA(1,1) with student's t 
% innovations (c,phi,theta,nu) and the respective
% time series y as input values.
% The function then calls the FUN_Likelihood_Students_t function to get the
% likelihood contributions in order to caluclate the Value o the Likelihood
% function. The Function returns the value of the Likelihood function
% Note: We changed the sign because many optimization algorithms are designed to minimize. 
% By changing the sign, maximizing the function becomes equivalent to minimizing its negative.


%Call FUN_Likelihood_Students_t function
Likelihood_contributions = FUN_Likelihood_Students_t(parameterset,y);

%Calculate Value of Likelihood Function
Value_LF = - sum(Likelihood_contributions);


end