function out=dynamico(X,k0,k1,omega0,approx)
% Computes the p-value of Onatski's dynamic test of H0: k=k0 vs. H1: k0<k<k1+1
% for the data in X (the cross-sectional dimension equals the number of rows,
% the time series dimension equals the number of columns).

% out=dynamico(X,k0,k1) compute p-value of the test using an automatically chosen 
% range of (relatively low) approximating frequences. No theoretical justification 
% for such a choice exists.

% out=dynamico(X,k0,k1,omega0) computes the p-value of the test for frequency omega0 
% If approx is not specified as an imput, the range of approximating frequences is 
% chosen around omega0 automatically. No theoretical justification for such a choice 
% exists.

% out=dynamico(X,k0,k1,omega0,approx) computes the p-value of the test using
% the approximating frequences of the form 2*pi*s/T with integers s specified in
% the column vector approx.
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
if k1-k0>18
    disp('Message: Asymptotic distribution of the test has not been tabulated for such a large k1')
    disp('Consider decreasing k1')
    return
end
if nargin<3
    disp('Input error: not enough arguments specified')
    return
elseif nargin==3
    omega0=0;
    approx=(4:min(floor(3.3*sqrt(T))+3,floor(T/2-0.1)))';
elseif nargin==4
    m0=floor(T*omega0/(2*pi));
    m1=min(m0,floor(3.3*sqrt(T)));
    approx=((m0-m1+4):min((m0-m1+3+floor(3.3*sqrt(T))),floor(T/2-0.1)))';
else
    [m,w]=size(approx);
    if m<w
        disp('Input error: approx must be a column vector')
        return
    end
    if sum(floor(approx)<approx)>0
        disp('Input error: approx must be a vector of integers')
        return
    end
    if m>T/2
        disp('Input error: The number of the approximating frequences must be smaller than T/2')
        return
    end
    if m<25
        disp('Warning: The asymptotic approximation may be poor.')
        disp('Consider increasing the number of approximating frequences')
    end
    if m<5*(k1-k0)
        disp('Warning: The asymptotic approximation may be poor.')
        disp('Consider increasing the number of approximating frequences.')
        disp('Alternatively, consider decreasing k1.')
    end
    for i=1:m
        for j=i:m
            if i==j
                if approx(i,1)==0|approx(i,1)==T/2
                    disp('Input error: Approximating frequences cannot be 0 or pi')
                    return
                end
            else
                if (approx(i,1)+approx(j,1))/T==floor((approx(i,1)+approx(j,1))/T)
                    disp('Input error: sum of any two approximating frequences cannot equal 0 modulo 2pi')
                    return
                elseif (approx(i,1)-approx(j,1))/T==floor((approx(i,1)-approx(j,1))/T)
                    disp('Input error: a difference of any two distinct approximating frequences cannot equal 0 modulo 2pi')
                    return
                end
            end
        end
    end
end
% Compute the dft's of the data
Xf=X*exp(-sqrt(-1)*2*pi*(1:T)'*approx'/T);
% Compute the test statistic
opts.disp=0;
evalues=abs(eigs(Xf'*Xf,k1+2,'lr',opts));
R=max((evalues(k0+1:k1,1)-evalues(k0+2:k1+1,1))./(evalues(k0+2:k1+1,1)-evalues(k0+3:k1+2,1)));
% Compute the p-value
out=sum(CVGUE(:,k1-k0)>R)/1000;