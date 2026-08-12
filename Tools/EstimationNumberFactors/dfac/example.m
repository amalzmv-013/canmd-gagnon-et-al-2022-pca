% clear

T=200;
N=100;
q=3;
s=1;
p=1;   % the dynamic factors are AR(1) processes
r=q*(s+1);
A=[.6; .3; 0];
randn('state',99);
rand('state',8);
maxit=100;
method=1; rmax=10; pmax=4;
stat1=[];stat2=[];rhat=[];
addpath ~/matlab/factors;
% DGP
for it=1:maxit;
u=randn(T+1+p,q);
f=u;
for j=1:q;
    f(:,j)=filter([1 A(j)],1,u(:,j));
end;
e=randn(T+p+1,N);
lambda=randn(N,q);
F=f;
for i=1:s;
  lambda=[lambda randn(N,q)];
  F=[F lagn(f,i)];
end;  
X=F*lambda'+e;
X=trimr(X,p+1,0);
F=trimr(F,p+1,0);
f=trimr(f,p+1,0);
% end DGP

if method==1; xepsilon=[2 2;1 1;.5 .5];end;     %covariance matrix
if method==2; xepsilon=[1.25 1.25;1 1;1.5 1.5];end;  % correlation matrix
eps=xepsilon/min([N^(.4);T^(.4)]);

[rhat(it),chat,Fhat,evals_s]=ICP(X,rmax,2,2);     % determine number of   static factros
[ehat,Fhat,lamhat,ve2]=pc(standard(X),rhat); % get rhat factors   using L'L/N normalization
[e,beta]=est_e(T,pmax,0,Fhat); % estimate VAR in Fht

[stat1(it,:),stat2(it,:),evals_d]=dfactest(e,N,T,eps,method);
end; % end it
