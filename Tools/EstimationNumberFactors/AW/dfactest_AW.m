function [qhat]=dfactest_AW(X,rmax,qmax,pmax,smax,stand);
% program to compute the number of dynamic factors
% and to compute U0 using SW method.
% Given estimate of F, F_hat, do the following regression:
% X_{i,t} = Ai(L)*X_{i,t-1} + Bi(L)*F_hat_{t} + eps_{i,t}
% Get residuals eps_hat_{t}. Given that 
% eps_{t} = Lambda*R*U0_{t} + e_{t}, 
% we can estimate the number of dynamic factors by applying Bai-Ng (2002)
% test, qo_hat, and then estimate U0 as qo_hat PC of eps_hat_{t}

% INPUTS:
% X - data
% rmax - number of static factors
% qmax - max number of dynamic factors to be tested
% pmax - AR lag order for X
% smax - lag order for Fhat

% OUTPUTS:
% V - estimate of U0, up to a rotation
% q00 - estimate of number of dynamic factors

% Dalibor Stevanovic
% October 14, 2011

%X=param;[T,N]=size(X);rmax=8;qmax=6;pmax=1;smax=1;stand=1;
% Step 0. Estimate rmax static factors from X
if stand==0
    xx=X;
elseif stand==1
    xx = standard(X);
end
[T,N]=size(xx);
%[Fhat,L] = extract_1(xx,rmax);  
[eee,Fhat,L,ss]=pc(xx,rmax);


% Step 1. Regress X on lags of X and lags of Fhat


dif = abs(size(xx,1)-size(Fhat,1));
if dif>0
    xx = trimr(xx,dif,0);
end
eps = NaN(T-dif,N);
for i=1:N
    reg1 = [];
    for p=1:pmax;
        reg1=[reg1 lagn(xx(:,i),p)];
    end;
    for s=1:smax;
        reg1=[reg1 lagn(Fhat,s)];
    end;
    y=xx(:,i);
    %reg=trimr(reg,pmax,0);
    AR=reg1\y;
    eps(:,i)=y-reg1*AR;
    
%     % NW regression
%     res=nwest(y,reg1,0);
%     eps(:,i)= y-res.yhat;
    
    clear reg1 y
end

% Step 2. Estimate the number of dynamic factors by applying a Bai-Ng(2002)
% test to eps
if stand==0
    DEMEAN = 0;
elseif stand==1
    DEMEAN = 2;
end
%disp(sprintf('Demean %d',DEMEAN));
%disp('detrmining number of factors');
%disp(sprintf('T= %d N= %d',size(eps)));
for i=1:7;
  %disp([nbpiid(eps,qmax,i,DEMEAN)   nbplog(eps,qmax,i,DEMEAN)]);
  qhat(i,:) = [nbpiid(eps,qmax,i,DEMEAN)   nbplog(eps,qmax,i,DEMEAN)];
end;
%q00 = nbplog(eps,qmax,2,2);

% % Step 3. Compute U0 as q0 PC of eps
% if q00>=q0
%     if stand==0
%         [V,L] = extract_1(eps,q00);
%     elseif stand==1
%         [V,L] = extract_1(standard(eps),q00);
%     end
% else
%     if stand==0
%         [V,L] = extract_1(eps,q0);
%     elseif stand==1
%         [V,L] = extract_1(standard(eps),q0);
%     end
% end
% %V = trimr(V,2,0); % to have the same dimension as in est_U.m file


