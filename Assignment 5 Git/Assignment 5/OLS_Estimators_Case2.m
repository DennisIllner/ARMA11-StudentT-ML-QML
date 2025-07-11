function [OLS_Estimates] = OLS_Estimators_Case2(y)

% The following Function performs an OLS regression for Case 2 of the Dickey-Fuller Test.
% The following regression is performed for a given time series y as inputvalue: y_t = alpha + roh*y_t-1 + e_t.
% The function returns the OLS estimate roh, its corresponding std. error, the value of
% the Dickeyfuller Teststatistic T*(roh-1) and the t-Value (roh-1)/s.e(roh)
% both under the null that roh = 1


%Perform OLS

y_reg=y(2:length(y)); % Dependent Variable (1 Obersavtion is lost since we regress on the previous value of y_t)

x_reg =y(1:length(y)-1); % Create "explainable Vector"
vector_one = ones(length(x_reg),1); % Vector for constant
x = [vector_one, x_reg]; % Merge both togehter in a matrix

b = (x'*x)^-1 * (x'*y_reg); %Calculation of regression coefficient (Here b is a 2x2 Variance-Covariance-Matrix)


%First calculating residual variance s^2
T = size(y_reg,1); %Use length of y_reg (true T since we lose one observation by the shift in the regression)
s_square = 1/(T-2) * sum((y_reg-b(1)-b(2)*x(:,2)).^2);

%Calculating estimated Variance of roh
Var_Cov_b = s_square * (x'*x)^-1;

%Calculating std. error of roh_hat
std_err_roh_hat =sqrt(Var_Cov_b(2,2));

%Calculate t-value
t_value = (b(2)-1)/sqrt(Var_Cov_b(2,2));

%Computing test statistic T(roh_hat-1)  
Dickey_Fuller = T*(b(2)-1)   ;

%Save Results
OLS_Estimates = [b(2),std_err_roh_hat,Dickey_Fuller,t_value];

end