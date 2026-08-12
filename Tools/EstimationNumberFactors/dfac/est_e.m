function [e,B]=est_e(T,p1,p2,Fhat);

% estimates a VAR in Fhat with p1 lags and P2 leads

r=cols(Fhat);
stat=zeros(1,r);
reg=ones(rows(Fhat),1);
B=[];
for j=1:p1;
reg=[reg lagn(Fhat,j)];
end;
for j=1:p2;
reg=[reg lagn(Fhat,-j)];
end;
Fhat=trimr(Fhat,p1+1,p2);
reg=trimr(reg,p1+1,p2);
%Fhat=trimr(Fhat,p1,p2);
%reg=trimr(reg,p1,p2);

e=Fhat;
for i=1:r;
  beta=reg\Fhat(:,i);
  B=[B beta];
  e(:,i)=Fhat(:,i)-reg*beta;
end;


