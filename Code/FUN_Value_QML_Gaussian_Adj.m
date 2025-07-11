function [Value_LF] = FUN_Value_QML_Gaussian_Adj(parameterset,y)
% The FUN_Value_QML_Gaussian_Adj function calculates the value of the 
% likelihood function under the assumption of a Gaussian (normal) distribution. 
% It takes the parameter set (c, phi, theta, sigma^2) and the data y.
%
% The function calls the FUN_QML_Gaussian function to obtain the individual 
% likelihood contributions and uses them to compute the overall value 
% of the likelihood function. The function returns the likelihood value.
% Note: We changed the sign because many optimization algorithms are designed to minimize. 
% By changing the sign, maximizing the function becomes equivalent to minimizing its negative.

%Call FUN_QML_Gaussian function
Likelihood_contributions = FUN_QML_Gaussian(parameterset,y);

%Calculate Value of Likelihood Function
Value_LF = - sum(Likelihood_contributions);


end