
%This program compute the number of factors in Ahn and Horenstein, 2013,
%Eigenvalue Ratio Test for the Number of Factors, Econometrica 2013
% This version is written by MAO TAKONGMO and is a modifying Gauss version written by Ahn and Horenstein

function  [ER,GR]=Ahn_Horenstein(X,rmax)

N=size(X,1);
T=size(X,2);

if N<T;
factors=X'*X/(T*N); 
else;
factors=X*X'/(T*N);
end;
a=eig(factors);
d=sort(a,1, 'descend');

% Eigenvalue criteria for choosing the number of factors */

m=min(N,T);

er=zeros(rmax,1);

for w=1:rmax;
er(w,1)=(d(w))/d(w+1);
end;
[C,storecrit]=max(er);

ds=zeros(rmax+1,1);

for i=1:rmax+1;
ds(i)=sum(d(i:min(N,T)))/sum(d(i+1:min(N,T)));
end;

gr=zeros(rmax+1,1);

for w=1:rmax;
gr(w,1)=log(ds(w))/log(ds(w+1));
end;

[C,storecrit2]=max(gr);
ER=storecrit;
GR=storecrit2;
