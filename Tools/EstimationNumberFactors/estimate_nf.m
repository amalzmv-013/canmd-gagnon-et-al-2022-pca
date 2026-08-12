function est_nf = estimate_nf(param,rmax,qmax,k0,k1,DEMEAN)

% MAP test
[m1,m2] = map_test(corr(param));
map_nf(1,:) = [m1,m2]; 
% map_nf: (squared part corr), (^4 part corr)
fprintf(1,'# static factors by MAP tests (2007)')
map_nf
fprintf(1,'\n')

% Parallel test
ndatsets = 1000; percent = 95; info=0;
iter=0;
for ii=1:2
    kind=ii;
    for jj=1:2
        iter=iter+1;
        randtype = jj;
        [numfac,realeval,means,percentiles]=parallel_test(param,ndatsets,percent,kind,randtype,info);
        paral_nf(1,iter)=numfac;
    end
end       
% paral: (kind=1,rand=1),(kind=1,rand=2),(kind=2,rand=1),(kind=2,rand=2)
fprintf(1,'# static factors by Parallel test')
paral_nf
fprintf(1,'\n')


% Estimate number of static factors Bai-Ng (2002)
fprintf(1,'# static factors by Bai-Ng (2002)')
disp(sprintf('Demean %d',DEMEAN));
disp('detrmining number of factors');
disp(sprintf('T= %d N= %d',size(param)));
for i=1:8;
  disp([nbpiid(param,rmax,i,DEMEAN)   nbplog(param,rmax,i,DEMEAN)]);
end;
fprintf(1,'\n')

% Estimate number of dynamic factors Bai-Ng(2007)
%[rhat,chat,Fhat,evals_s]=ICP(param,rmax,7,1);     % determine number of   static factros
[ehat,Fhat,lamhat,ve2]=pc(standard(param),rmax); % get rhat factors   using L'L/N normalization
[e,beta]=est_e(size(param,1),qmax,0,Fhat); % estimate VAR in Fht
xepsilon=[2 2;1 1;.5 .5];     %covariance matrix
eps=xepsilon/min([size(param,2)^(.4);size(param,1)^(.4)]);
[stat1,stat2,evals_d]=dfactest(e,size(param,2),size(param,1),eps,1);
fprintf(1,'# dynamic factors by Bai-Ng (2007)')
[stat1]
[stat2]
fprintf(1,'\n')

% Estimate number of dynamic factors Amengual-Watson(2007)
fprintf(1,'# dynamic factors by Amengual-Watson (2007)')
[qhat]=dfactest_AW(param,rmax,qmax,1,1,1)
fprintf(1,'\n')

% Estimate number of factors Hallin and Liska (2007)
[nfact, v_nfact, cr] = numfactors(param,rmax,floor(size(param,2)/4),1,1);

%% Estimate number of factors ABC (2010)
ABC_crit(param,30,1,floor(size(param,2)/4));

% Estimate number of factors Onatski (2010)
fprintf(1,'# static factors by Onatski (2010)')
onatski_k = onatski(standard(param)',rmax)
fprintf(1,'\n')

% Estimate number of factors Onatski (2009)
% static: H0: k=k0 vs. H1: k0<k<k1+1
for k0i = 1:k0
    for k1i = 1:k1
        onatski_k_09(k0i,k1i) = statico(standard(param)',k0i,k0i+k1i);
    end
end
fprintf(1,'# static factors by Onatski (2009) \n')
fprintf(1,'H0: k=k0 vs. H1: k0<k<k1+1 \n')
fprintf(1,'k0 lines, k1 columns \n')
onatski_k_09
fprintf(1,'\n')
% dynamic
for k0i = 1:k0
    for k1i = 1:k1
        onatski_q_09(k0i,k1i) = dynamico(standard(param)',k0i,k0i+k1i);
    end
end
fprintf(1,'# dynamic factors by Onatski (2009) \n')
fprintf(1,'H0: k=k0 vs. H1: k0<k<k1+1 \n')
fprintf(1,'k0 lines, k1 columns \n')
onatski_q_09
