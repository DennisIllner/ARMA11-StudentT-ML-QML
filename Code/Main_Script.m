%%%%%%%%%%%%%%%%%%%%%%%%
%%Mandatory Assignment%%
%%%%%%%%%%%%%%%%%%%%%%%%
%This is the main File%

clc
clear
close all
rng(415) %set seed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Task 1) Simulating ARMA(1,1)

%1.1 Create Function


%1.2 %Call Function

%Specifiy Parameters
c=2;
phi = 0.95;
theta = 0.25;
nu = 4;
y_0 = c/(1-phi); %Mean of ARMA(1,1) 
T=800;

y_800 = FUN_ARMA11(T,c,phi,theta,nu,y_0);

y_800 = y_800(51:T) ; %Drop first fifty obersvations as burn-in phase avoiding dependencies on the first obs.


%1.3 Plot one Realization of the Process
y_mean = mean(y_800,1);
t = (1:length(y_800))'; %Create discrete x-axis for time t


%Figure for Analyzes
figure
handle = plot(t, y_800, t, y_mean); %crerate handle for formatting purposes
handle2 = yline(y_mean, 'r--', 'LineWidth', 2);
legend('$y_t$', 'Location', 'NorthWest', 'Interpreter', 'latex'); %create legend
xlabel('$\mathbf{t}$', 'FontSize', 14, 'Interpreter', 'latex', 'FontWeight', 'bold'); %describe axis
ylabel('$\mathbf{y}$', 'Interpreter', 'latex', 'FontWeight', 'bold'); %describe axis
xlim([0, 750]) %Scale x-axis
title('$\mathbf{Realization\ of\ an\ ARMA(1,1)\ with\ c=2,\ \phi =0.95,\ \theta =0.25,\ \nu =4,\ y_0=\frac{c}{1-\phi}}$', 'Interpreter', 'latex', 'FontWeight', 'bold');  
set(handle(1), 'Color', 'blue', 'LineWidth', 2)  %style-Formatting

%Figure for paper
figure
handle = plot(t, y_800, t, y_mean); % Create handle for formatting purposes
handle2 = yline(y_mean, 'r--', 'LineWidth', 2); % Dashed line for the mean
legend([handle(1), handle2], '$\mathbf{y_t}$', '\textbf{mean}', 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('$\mathbf{t}$', 'FontSize', 20, 'Interpreter', 'latex', 'FontWeight', 'bold'); 
ylabel('$\mathbf{y}$', 'FontSize', 20, 'Interpreter', 'latex', 'FontWeight', 'bold'); 
xlim([0, 750]); % Scale x-axis
set(handle(1), 'Color', 'blue', 'LineWidth', 2); % Style formatting
set(gca, 'FontSize', 16, 'FontWeight', 'bold'); % Thicker font for axis 
set(gca, 'Box', 'on', 'LineWidth', 2); % Box around plot



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Task2 Unit Root Tests

%2.1 %Decide which Case to choose for the Dickey Fuller Test

disp('Nr.2.1');
disp('Anaylising the time series from Task 1.3 we can see a mean reverting behaviour around approx. 39.');
disp('We can not detect a trend in the series.');
disp('Case 2 with y_t = a + roh*y_t-1 + et seems to be the best case for the dickey-fuller test.');
disp(' ')

%2.2 / 2.3 / 2.4

%Using the Function for the unit root test of Case 2 to compute OLS
%estimator, OLS standard errors and test-statistic.
results =OLS_Estimators_Case2(y_800);

%Extract relevant estimates
OLS_estimator = results(1,1)  ;        %OLS estimator
OLS_standard_error = results(1,2);     %OLS standard error
DF_Teststatistic = results(1,4) ;      %t-statistic

%Extract Critical Value T_crit = 1- alpha where alpha 5% (Hamilton)
t_crit = -2.86;

%Interpretation
disp('Nr. 2.2/2.3/2.4')
disp(['OLS estimator: ' num2str(OLS_estimator)])
disp(['OLS standard error: ' num2str(OLS_standard_error)])
disp(['DF_Teststatistic: ' num2str(DF_Teststatistic)])
disp(['Critical Value: ' num2str(t_crit)])
disp(' ')
disp('We can reject the null hypothesis H0: p=1 (presence of a unit root) at the 5% significance level using the t-statistic.');
disp(['The t-statistic is ' num2str(DF_Teststatistic) ', which is smaller than the critical value t_crit = ' num2str(t_crit) '.']);
disp('Therefore, the test statistic falls in the rejection region.');
disp('As a result, we find sufficient evidence to reject the null hypothesis of a unit root at the 5% significance level.');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Nr.3 % Conditional Likelihood

%3.1 Create function calculating log-likelihood contributions
%Function is called FUN_Likelihood_Students_t


%3.2 Create another Function which calls FUN_Likelihood_Students_t and calculates
%the value of the function for a given paramterset and time series y


%a)
%Specify Parameters
c     = 2;
phi   = 0.95;
theta = 0.25;
nu    = 4;


%Save parameters in vector
parameters =[c,phi,theta,nu];


%Call FUN_Value_LF_Students_t function to get the Values of the LF
Value_a = FUN_Value_LF_Students_t(parameters,y_800);


%b)
%Specify Parameters
c     = 1.5;
phi   = 0.75;
theta = 0.5;
nu    = 6;


%Save parameters in vector
parameters_wrong =[c,phi,theta,nu];


%Call FUN_Value_LF_Students_t function to get the Values of the LF
Value_b = FUN_Value_LF_Students_t(parameters_wrong,y_800);


%Interpretation
disp(' ')
disp('Nr.3.2)')
disp(['Log-Likelihood Contributions for a): ' num2str(Value_a)]);
disp(['Log-Likelihood Contributions for b): ' num2str(Value_b)]);
disp('The value of the likelihood function for a) is higher than that for b).');
disp('The value of the likelihood function indicates how well the parameters fit the data of the ARMA(1,1) process.');
disp('A higher likelihood value suggests that the parameters provide a better fit to the data.');
disp('Therefore, we conclude that the parameters in a) describe the ARMA(1,1) data better than those in b).');
disp('We know that the parameters in a) are the true parameters of the ARMA(1,1) process.');
disp('Therefore, these parameters yield a high value of the likelihood function,');
disp('and thus the likelihood value for a) should be greater than the value for b) when using paramters which deviate from the true parameters.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Task 4

%4.1 Include th CML toolbox
addpath('CML');

%4.2 % Parameter estimation using CML toolbox


%1.First adjust Function from 3.2 such that minimizing is equal
%to maximizing the function (change sign)
%We then get the new function called: FUN_Value_LF_Students_t_Adj


%Starting Values
x0 = [1.5,0.75,0.5,5]';


%Set optimization options
options = optimset('Display','off','TolX',10 ^(-40),'TolFun',10 ^(-40),'MaxIter',10^10 , 'MaxFunEvals', 100000); 


%Call Function
%Include Inputvalues: FUN_Value_LF_Students_t_Adj, FUN_Likelihood_Students_t, data y,
%starting values x0, fminsearch (1), Hessen-based estimation (1) and options
[x,f,g,cov,retcode]= CML(@FUN_Value_LF_Students_t_Adj, @FUN_Likelihood_Students_t, y_800, x0, 1, 1, options);



%4.3 Computing Standard Errors of Paramterestimates and their 95%
%Confidence Intervalls


%Compute Standard Errors via the Var-Cov-Matrix of the estimates. Extract
%relevant variances and square them.
variances = diag(cov);
std_errors = sqrt(variances);


%Compute Confidence Intervalls
% alpha = 5% Significance level. Quantile of the std. normal distribution 
% at 1-( alpha/2)
quantile = norminv(0.975, 0, 1);


%Create Confidence Intervals [theta_hat -/+ s.e(theta_hat)*quantile]
Lower_Bounds = x-std_errors*quantile;
Upper_Bounds = x+std_errors*quantile;
Confidence_Intervals095 = [Lower_Bounds,Upper_Bounds];


% Create a table to illustrate results
parameter_names = {'c', 'phi', 'theta', 'nu'}';
results_table = table(parameter_names, x, std_errors, Lower_Bounds, Upper_Bounds,'VariableNames', {'Parameter', 'Estimate', 'Std_Error', 'Lower_Bound', 'Upper_Bound'});


% Display the results and interpretation
disp(' ')
disp('Nr.4.3')
disp('Maximum Likelihood estimates, it`s standard errors, and 95% confidence intervals for T=750:');
disp(results_table);
disp(' ')
disp('General Interpretation:')
disp('With 95% of the cases, the true value of the respective parameters lies within the confidence bounds.');
disp('The confidence interval is constructed such that it contains the true value in 95% of repeated samples.');
disp('For all values inside the interval constructed with the parameter estimates, I cannot reject the null hypothesis H0: theta_0 = theta_bar,');
disp('for chosen values of theta_bar that lie within this interval.');


%4.4 % Perform Hypothesis Test
disp(' ')
disp('Nr.4.4')
disp('Performing Hypothesis test for parameters:');
disp('1: Define Hypithesis: H0: phi = 0.8   HA: phi ~= 0.8')
disp('2: Siginificance Level is alpha = 5%');
disp(['3: Determine non-rejection area (two sided): [' num2str(-norminv(0.975, 0, 1)) ',' num2str(norminv(0.975, 0, 1)) ']']);
disp(['4: Caluculate t-value: ' num2str( (x(2)-0.8)/std_errors(2) )]);
disp(['5: Interpretation:  Since ' num2str((x(2)-0.8)/std_errors(2)) ' > ' num2str(norminv(0.975, 0, 1)) ', we  can reject' ]);
disp('the null hypothesis that phi=0.8 with a given significance level of 5%.');
disp('This also confirms the statement that the confidence interval provides all values for theta_bar for which we cannot reject H0.');
disp('Since theta_bar = 0.8 lies outside the calculated interval above, we can reject H0 at the 5% significance level.');


%4.5 % Compute p_value
t_val = (x(2)-0.8)/std_errors(2);
p_value = 2 * (1 - normcdf(t_val, 0, 1)) ;

disp(' ')
disp('Nr.4.5')
disp(['Interpretation: The p-value of ' num2str(p_value) ' indicates that we can reject the null hypothesis at any conventional significance level.']) 
disp('This means we find strong evidence that the true value is highly significantly different from 0.8.');



%4.6 Performing Task 4.2-4.5 using T=50000


%Specifiy Parameters and generate data of an ARMA(1,1) using T=50000
c=2;
phi = 0.95;
theta = 0.25;
nu = 4;
y_0 = c/(1-phi); %Mean of ARMA(1,1) 
T=50000;


%Generate Data
y_50k = FUN_ARMA11(T,c,phi,theta,nu,y_0); 


%Drop first fifty obersvation as burn-in phase avoiding dependencies on the first obs.
y_50k = y_50k(51:T) ; 


%Call CML Function with respective starting values and options
x0 = [1.5,0.75,0.5,5]';

options = optimset('Display','off','TolX',10 ^(-40),'TolFun',10 ^(-40),'MaxIter',10^10 , 'MaxFunEvals', 100000); 

[x,f,g,cov,retcode]= CML(@FUN_Value_LF_Students_t_Adj, @FUN_Likelihood_Students_t, y_50k, x0, 1, 1, options);


%Extract variances and std. errors of parameter estimates
variances = diag(cov);
std_errors = sqrt(variances);


%Compute quantile to construct confidence intervals
quantile = norminv(0.975, 0, 1);


%Create Confidence Intervals [theta_hat -/+ s.e(theta_hat)*quantile]
Lower_Bounds = x-std_errors*quantile;
Upper_Bounds = x+std_errors*quantile;
Confidence_Intervals095 = [Lower_Bounds,Upper_Bounds];


% Create a table to illustrate results
parameter_names = {'c', 'phi', 'theta', 'nu'}';
results_table = table(parameter_names, x, std_errors, Lower_Bounds, Upper_Bounds,  ...
    'VariableNames', {'Parameter', 'Estimate', 'Std_Error', 'Lower_Bound', 'Upper_Bound'});


% Display the results and interpretation
disp(' ')
disp('Nr.4.6')
disp('Maximum Likelihood estimates, it`s standard errors, and 95% confidence intervals for T=49950:');
disp(results_table);


% Perform Hypothesis Test
disp(' ')
disp('Performing Hypothesis test for parameters:');
disp('1: Define Hypithesis: H0: phi = 0.8   HA: phi ~= 0.8')
disp('2: Siginificance Level is alpha = 5%');
disp(['3: Determine non-rejection area (two sided): [' num2str(-norminv(0.975, 0, 1)) ',' num2str(norminv(0.975, 0, 1)) ']']);
disp(['4: Calculate t-value: ' num2str( (x(2)-0.8)/std_errors(2) )]);
disp(['5: Interpretation:  Since ' num2str((x(2)-0.8)/std_errors(2)) ' > ' num2str(norminv(0.975, 0, 1)) ', we  can reject' ]);
disp('the null hypothesis that phi=0.8 with a given significance level of 5%.');
disp('This also confirms the statement that the confidence interval provides all values for theta_bar for which we cannot reject H0.');
disp('Since theta_bar = 0.8 lies outside the calculated interval above, we can reject H0.');


% Compute p_value
t_val = (x(2)-0.8)/std_errors(2);
p_value = 2 * (1 - normcdf(t_val, 0, 1)) ;
disp(' ')
disp(['Interpretation p-value: The p-value of ' num2str(p_value) ' indicates that we can reject the null hypothesis at any conventional significance level.']) 
disp('This means we find strong evidence that the true value is highly significantly different from 0.8.');


%Interpretation 4.6
disp(' ')
disp('Interpretation under consideration of change in T and covariance matrix:')
disp('With an increased sample size T, we observe that the parameter estimates converge closer to the true values.')
disp('This reflects the consistency of the estimators, which ensures that as T grows, the estimates converge in probability to the true parameters.') 
disp('The variances, and thus the standard errors, are smaller, resulting in narrower confidence intervals.') 
disp('Regarding hypothesis testing, we obtain much higher t-values, making it more likely to reject H0. Consequently, the p-values become smaller.') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Nr.5 Quasi Maximum Likelihood

%5.1 Write Functions to calculate Quasi Maximum Likelihood Contributions
%and Value of LF (FUN_QML_Gaussian and FUN_Value_QML_Gaussian)




%Nr.5.2
disp(' ')
disp('Nr.5.2')
disp('Quasi-Maximum-Likelihood (QML) Estimation:');
disp('QML estimation applies when we are uncertain about the underlying distribution,');
disp('and assume a potentially incorrect distribution for the maximization of the Likelihood Function.');
disp('However when using QML we rely on the Law of Large Numbers (LLN),');
disp('which holds if our data comes from a stationary and ergodic process. We further assume that the limit function is');
disp('maximized at the pseudo-true parameters. Under these conditions, we obtain consistent an asymptotically normal distributed estimators.');
disp('An incorrect distribution assumption may lead us not to estimate the "true" parameters.');
disp('However, the estimators are CAN (consistent and asymptotically normal).');
disp('Interpretation under an incorrect distribution assumption:');
disp('1.QML minimizes the Kullback-Leibler Criteria.');
disp('2.If the first moment of the conditional mean is correctly specified,');
disp('but the conditional density is wrong, the parameters defining the conditional mean are consistent');



%Nr.5.3
%T=800

%Call CML Function with respective starting values and options
x0 = [1.5,0.75,0.5,1]';

options = optimset('Display','off','TolX',10 ^(-40),'TolFun',10 ^(-40),'MaxIter',10^10 , 'MaxFunEvals', 100000); 

[x_QML_800,f,g,cov_QML_800,retcode]= CML(@FUN_Value_QML_Gaussian_Adj, @FUN_QML_Gaussian, y_800, x0, 1, 3, options); %Adjust COVPAR option to QML-based


%Nr.5.4 Computing Standard Errors of parameter estimates and their 95%
%Confidence Intervalls

%Extract variances and std. errors of parameter estimates
variances_QML_800 = diag(cov_QML_800);
std_errors_QML_800 = sqrt(variances_QML_800);

%Compute 95% Confidence Intervales
quantile = norminv(0.975, 0, 1);


%Create Confidence Intervals [theta_hat -/+ s.e(theta_hat)*quantile]
Lower_Bounds_QML_800 = x_QML_800-std_errors_QML_800*quantile;
Upper_Bounds_QML_800 = x_QML_800+std_errors_QML_800*quantile;
Confidence_Intervals095_QML_800 = [Lower_Bounds_QML_800,Upper_Bounds_QML_800];

% Create a table to illustrate results
parameter_names = {'c', 'phi', 'theta', 'sigma^2'}';
results_table = table(parameter_names, x_QML_800, std_errors_QML_800, Lower_Bounds_QML_800, Upper_Bounds_QML_800,'VariableNames', {'Parameter', 'Estimate', 'Std_Error', 'Lower_Bound', 'Upper_Bound'});


% Display the results and interpretation
disp(' ')
disp('Nr.5.4')
disp('Quasi Maximum Likelihood estimates, it`s standard errors, and 95% confidence intervals for T=750:');
disp(results_table);


%%%%%%%%
%T50000

%Call CML Function with respective starting values and options
x0 = [1.5,0.75,0.5,1]';

options = optimset('Display','off','TolX',10 ^(-40),'TolFun',10 ^(-40),'MaxIter',10^10 , 'MaxFunEvals', 100000); 

[x_QML_50k,f,g,cov_QML_50k,retcode]= CML(@FUN_Value_QML_Gaussian_Adj, @FUN_QML_Gaussian, y_50k, x0, 1, 3, options); %Adjust COVPAR option to QML-based


%Nr.5.4 Computing Standard Errors of Paramterestimates and their 95%
%Confidence Intervalls

%Extract variances and std. errors of parameter estimates
variances_QML_50k = diag(cov_QML_50k);
std_errors_QML_50k = sqrt(variances_QML_50k);

%Compute 95% Confidence Intervales
quantile = norminv(0.975, 0, 1);


%Create Confidence Intervals [theta_hat -/+ s.e(theta_hat)*quantile]
Lower_Bounds_QML_50k = x_QML_50k-std_errors_QML_50k*quantile;
Upper_Bounds_QML_50k = x_QML_50k+std_errors_QML_50k*quantile;
Confidence_Intervals095_QML_50k = [Lower_Bounds_QML_50k,Upper_Bounds_QML_50k];

% Create a table to illustrate results
parameter_names = {'c', 'phi', 'theta', 'sigma^2'}';
results_table = table(parameter_names, x_QML_50k, std_errors_QML_50k, Lower_Bounds_QML_50k, Upper_Bounds_QML_50k,'VariableNames', {'Parameter', 'Estimate', 'Std_Error', 'Lower_Bound', 'Upper_Bound'});


% Display the results and interpretation
disp(' ')
disp('Quasi Maximum Likelihood estimates, it`s standard errors, and 95% confidence intervals for T=49950:');
disp(results_table);

disp(' ')
disp('For Larger T the QML estimates converge to the true parameters except the')
disp('parameter for sigma^2 which is not defined in the conditional mean')





%Nr.6

%Nr.6.1

%Set parameters
c=2;
phi = 0.95;
theta = 0.25;
nu = 4;
y_0 = c/(1-phi); 
T=800;
n = 2500; %numer of ensembles
burn_in = 50; %burn-in phase


%Initialize Matrix for ensembles
ensembles = zeros((T-burn_in),n);


%Initialize 4 matrices to save parameter estimates and std_errors for ML and QML 
ML_estimates = zeros(n,4);
ML_estimates_std_errors = zeros(n,4);
QML_estimates = zeros(n,4);
QML_estimates_std_errors = zeros(n,4);

%Parallel Computing to increase computational Power
%Parallelize to optimize simulation
%numCores = feature('numCores'); % Get the number of available cores
%parpool(numCores - 1); % Open the parallel pool with all cores except one 

%Simulation
for i=1:size(ML_estimates,1)

%1.Generate Date
y = FUN_ARMA11(T,c,phi,theta,nu,y_0);

y = y((burn_in+1):T) ; %Drop first fifty obersvations as burn-in phase avoiding dependencies on the first obs.

ensembles(:,i)= y;


%2. Perform ML estimation 

%Call CML Function with respective starting values and options
x0 = [1.5,0.75,0.5,5]';
options = optimset('Display','off','TolX',10 ^(-40),'TolFun',10 ^(-40),'MaxIter',10^10 , 'MaxFunEvals', 100000); 
[x,f,g,cov,retcode]= CML(@FUN_Value_LF_Students_t_Adj, @FUN_Likelihood_Students_t, ensembles(:,i), x0, 1, 1, options);

%Save estimates and std_errors
ML_estimates(i,:)=x';
std_errors_ML = sqrt(diag(cov));
ML_estimates_std_errors(i,:) = std_errors_ML;


%3. Perform QML estimation

%Call CML Function with respective starting values and options
x0 = [1.5,0.75,0.5,1]';
options = optimset('Display','off','TolX',10 ^(-40),'TolFun',10 ^(-40),'MaxIter',10^10 , 'MaxFunEvals', 100000); 
[x_QML,f,g,cov_QML,retcode]= CML(@FUN_Value_QML_Gaussian_Adj, @FUN_QML_Gaussian, ensembles(:,i), x0, 1, 3, options); %Adjust COVPAR option to QML-based

%Save estimates and std_errors
QML_estimates(i,:)=x_QML';
std_errors_QML = sqrt(diag(cov_QML));
QML_estimates_std_errors(i,:) = std_errors_QML;
 
end

%delete(gcp('nocreate'))   % Delete old pool


%6.2 Kernel Densities of parameters for ML and QML

%Compute densities of respective parameters 
[density_c_ML,x_c_ML]=ksdensity(ML_estimates(:,1));
[density_c_QML,x_c_QML]=ksdensity(QML_estimates(:,1));

[density_phi_ML,x_phi_ML]=ksdensity(ML_estimates(:,2));
[density_phi_QML,x_phi_QML]=ksdensity(QML_estimates(:,2));

[density_theta_ML,x_theta_ML]=ksdensity(ML_estimates(:,3));
[density_theta_QML,x_theta_QML]=ksdensity(QML_estimates(:,3));



%Figure for analyzes
%Create Subplots for each parameter estimate
figure;
subplot(3,1,1)
handle = plot(x_c_ML, density_c_ML, x_c_QML, density_c_QML); 
legend('$\hat{c}$ ML ','$\hat{c}$ QML',  'Location', 'NorthWest', 'Interpreter', 'latex');
title('Distribution of $\hat{c}$ for ML and QML ','Interpreter', 'latex')
ylabel('density', 'FontSize', 20, 'Interpreter', 'latex', 'FontWeight', 'bold');
set(handle(1), 'Color', 'red', 'LineWidth', 1);
set(handle(2), 'Color', 'blue', 'LineWidth', 1);


subplot(3,1,2)
handle = plot(x_phi_ML, density_phi_ML, x_phi_QML, density_phi_QML); 
legend('$\hat{\phi}$ ML ','$\hat{\phi}$ QML',  'Location', 'NorthWest', 'Interpreter', 'latex');
title('Distribution of $\hat{\phi}$ for ML and QML ','Interpreter', 'latex')
ylabel('density', 'FontSize', 20, 'Interpreter', 'latex', 'FontWeight', 'bold');
set(handle(1), 'Color', 'red', 'LineWidth', 1);
set(handle(2), 'Color', 'blue', 'LineWidth', 1);


subplot(3,1,3)
handle = plot(x_theta_ML, density_theta_ML, x_theta_QML, density_theta_QML); 
legend('$\hat{\theta}$ ML ','$\hat{\theta}$ QML',  'Location', 'NorthWest', 'Interpreter', 'latex');
title('Distribution of $\hat{\theta}$ for ML and QML ','Interpreter', 'latex')
ylabel('density', 'FontSize', 20, 'Interpreter', 'latex', 'FontWeight', 'bold');
set(handle(1), 'Color', 'red', 'LineWidth', 1);
set(handle(2), 'Color', 'blue', 'LineWidth', 1);


%Figure for paper
%Create Subplots for each parameter estimate
figure;
subplot(3,1,1)
handle = plot(x_c_ML, density_c_ML, x_c_QML, density_c_QML); 
legend('$\hat{c}$ ML ','$\hat{c}$ QML',  'Location', 'NorthWest', 'Interpreter', 'latex');
ylabel('density', 'FontSize', 20, 'Interpreter', 'latex', 'FontWeight', 'bold');
set(handle(1), 'Color', 'red', 'LineWidth', 2);
set(handle(2), 'Color', 'blue', 'LineWidth', 2);
set(gca, 'FontSize', 16, 'FontWeight', 'bold'); % Thicker font for axis tick numbers and axis labels
set(gca, 'Box', 'on', 'LineWidth', 2); % Box around plot




subplot(3,1,2)
handle = plot(x_phi_ML, density_phi_ML, x_phi_QML, density_phi_QML); 
legend('$\hat{\phi}$ ML ','$\hat{\phi}$ QML',  'Location', 'NorthWest', 'Interpreter', 'latex');
ylabel('density', 'FontSize', 20, 'Interpreter', 'latex', 'FontWeight', 'bold');
set(handle(1), 'Color', 'red', 'LineWidth', 2);
set(handle(2), 'Color', 'blue', 'LineWidth', 2);
set(gca, 'FontSize', 16, 'FontWeight', 'bold'); 
set(gca, 'Box', 'on', 'LineWidth', 2); 


subplot(3,1,3)
handle = plot(x_theta_ML, density_theta_ML, x_theta_QML, density_theta_QML); 
legend('$\hat{\theta}$ ML ','$\hat{\theta}$ QML',  'Location', 'NorthWest', 'Interpreter', 'latex');
ylabel('density', 'FontSize', 20, 'Interpreter', 'latex', 'FontWeight', 'bold');
set(handle(1), 'Color', 'red', 'LineWidth', 2);
set(handle(2), 'Color', 'blue', 'LineWidth', 2);
set(gca, 'FontSize', 16, 'FontWeight', 'bold'); 
set(gca, 'Box', 'on', 'LineWidth', 2); 


%6.3
disp(' ')
disp('Nr.6.3')
disp('Kernel densitys for the parameter estimates of ML und QML are centered around the true paramter values,')
disp('reflecting the consistency property. Moreover, the kernel densities of the ML estimators are narrower,')
disp('indicating lower variances. This observation aligns with the efficiency property of the ML estimates.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Nr.7 Stochastic Confidence Intervals


%Maximum Likelihood

%Determine quantile
quantile = norminv(0.975, 0, 1);

% Initialise matrix for intervalls
Confidence_Intervals095 = zeros(size(ML_estimates, 1), size(ML_estimates, 2) * 2); 


% Calculate lower and upper bounds for respective parameters
for i = 1:size(ML_estimates, 2)
    Lower_Bounds = ML_estimates(:, i) - ML_estimates_std_errors(:, i) * quantile;
    Upper_Bounds = ML_estimates(:, i) + ML_estimates_std_errors(:, i) * quantile;
    
    Confidence_Intervals095(:, 2*i-1) = Lower_Bounds;
    Confidence_Intervals095(:, 2*i) = Upper_Bounds;

end


% Initialise count matrix
count_intervals = zeros(n,size(ML_estimates, 2) ); 

% Iterate over respective intervalls
for i = 1:4
   
    %Take respective interval for parameter we look at
    lower_bound = Confidence_Intervals095(:, 2*i-1); 
    upper_bound = Confidence_Intervals095(:, 2*i);     

    % Check if parametervalue from 1.1 is inside the intervall and assign
    % value 1 or 0
    count_intervals(:, i) = count_intervals(:, i) + ((parameters(i) >= lower_bound) & (parameters(i) <= upper_bound));
end

sum_inside = sum(count_intervals,1);

fraction = sum_inside/n;


% Create a table to illustrate results
parameter_names = {'c','phi','theta','nu'}';
results_table = table(parameter_names, sum_inside', fraction','VariableNames',{'Parameter','Number Inside','Fraction'});


disp(' ')
disp('Nr.7.1')
disp('Stochastic Confidence Intervals ML ')
disp(results_table);
disp('The frequencies of the true parameter values falling within the 95% confidence intervals are closely to 95%. This results from the consistency and asymptotic normality property of the ML estimates')

%%%%Not asked%%%%      
% Quasi Maximum Likelihood

% Initialise matrix for intervall
Confidence_Intervals095_QML = zeros(size(QML_estimates, 1), size(QML_estimates, 2) * 2); 


% Calculate lower and upper bounds for respective parameters
for i = 1:size(QML_estimates, 2)
    Lower_Bounds = QML_estimates(:, i) - QML_estimates_std_errors(:, i) * quantile;
    Upper_Bounds = QML_estimates(:, i) + QML_estimates_std_errors(:, i) * quantile;
    
    Confidence_Intervals095_QML(:, 2*i-1) = Lower_Bounds;
    Confidence_Intervals095_QML(:, 2*i) = Upper_Bounds;

end


% Initialise count matrix
count_intervals = zeros(n,size(QML_estimates, 2) ); 

% Iterate over respective intervalls
for i = 1:4
   
    %Take respective intervall for parameter we look at
    lower_bound = Confidence_Intervals095_QML(:, 2*i-1); 
    upper_bound = Confidence_Intervals095_QML(:, 2*i);     

    % Check if parametervalue from 1.1 is inside the intervall and assign
    % value 1 or 0
    count_intervals(:, i) = count_intervals(:, i) + ((parameters(i) >= lower_bound) & (parameters(i) <= upper_bound));
end

sum_inside = sum(count_intervals,1);

fraction = sum_inside/n;


% Create a table to illustrate results
parameter_names = {'c','phi','theta','sigma^2'}';
results_table_QML = table(parameter_names, sum_inside', fraction','VariableNames',{'Parameter','Number Inside','Fraction'});


disp(' ')
disp('Nr.7.1')
disp('Stochastic Confidence Intervals QML ')
disp(results_table_QML);


