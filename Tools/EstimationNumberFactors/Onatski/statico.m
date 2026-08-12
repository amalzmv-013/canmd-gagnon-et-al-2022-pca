function out=statico(X,k0,k1)
% Computes the p-value of Onatski's static test of H0: k=k0 vs. H1: k0<k<k1+1
% for the data in X (the cross-sectional dimension equals the number of rows,
% the time series dimension equals the number of columns).

load CVGUE.dat
[n,T]=size(X);
if k1-k0>7
    disp('Warning: Power of the test may be low')
    disp('Consider decreasing k1.')
end
if n<30
    disp('Warning: The asymptotic approximation may be poor.')
    disp('The cross-sectional dimension of the data is very low')
end
if T<50
    disp('Warning: The asymptotic approximation may be poor.')
    disp('The time-series dimension of the data is very low')
end
if k1-k0>18
    disp('Message: Asymptotic distribution of the test has not been tabulated for such a large k1')
    disp('Consider decreasing k1')
    return
end
if nargin<3
    disp('Input error: not enough arguments specified')
    return
end
% Compute the transformation of the data
T1=floor(T/2);
Xf=X(:,1:T1)+sqrt(-1)*X(:,(T1+1):(2*T1));
% Compute the test statistic
opts.disp=0;
%evalues=abs(eigs(Xf'*Xf,k1+2,'lr',opts)); % Matlab 2014 eigs function returns eigenvalues in ascending order!
evalues=abs(flipud(eigs(Xf'*Xf,k1+2,'lr',opts)));
R=max((evalues(k0+1:k1,1)-evalues(k0+2:k1+1,1))./(evalues(k0+2:k1+1,1)-evalues(k0+3:k1+2,1)));
% Compute the p-value
out=sum(CVGUE(:,k1-k0)>R)/1000;