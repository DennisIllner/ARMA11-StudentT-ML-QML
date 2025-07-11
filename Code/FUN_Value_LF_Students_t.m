function [Value_LF] = FUN_Value_LF_Students_t(parameterset,y)
% The FUN_Value_LF_Students_t function calculates the value of the
% Likelihood function for the student's t distribution
% and takes the parameterset of an ARMA(1,1) with student's t 
% innovations (c,phi,theta,nu) and the respective
% time series y as input values.
% The function then calls the FUN_Likelihood_Students_t function to get the
% likelihood contributions in order to caluclate the Value o the Likelihood
% function. The Function returns the value of the Likelihoodfunction


%Call FUN_Likelihood_Student's t function
Likelihood_contributions = FUN_Likelihood_Students_t(parameterset,y);

%Calculate Value of Likelihood Function
Value_LF = sum(Likelihood_contributions);


end