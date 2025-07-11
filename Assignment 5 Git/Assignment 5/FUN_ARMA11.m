function [y] = FUN_ARMA11(T,c,phi,theta,nu,y_0)
%The FUN_ARMA11 Function generates Data following an ARMA(1,1) process. The
%input parameters are
% T - sample size
% c - constant
% phi - Ar(1) Coefficient
% theta - MA(1) Coefficient
% nu - Degress of freedom t-distribution (nu)
% y_0 - Starting Value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Generate t-distributed innovations e with degrees of freedom nu
%T+1 to ensure we have one extra value for e(0)
e = trnd(nu,T+1,1) ;

%Intitalize output vector
y = zeros(T,1) ;

%Create first value 
y(1) = c + phi*y_0 + e(2) + theta*e(1) ;

%Using Loop to iterate forward and generate the process
for i = 2:T
    y(i) = c + phi*y(i-1) + e(i+1) + theta*e(i) ;

end


end