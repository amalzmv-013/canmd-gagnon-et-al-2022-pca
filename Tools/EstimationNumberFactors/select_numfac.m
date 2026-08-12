function [estim_numfac] = select_numfac(X,rmax)

%X = QC; rmax = 20;
[T,N] = size(X);

X = standard(X);

    % Bai-Ng (2002)
    for i=1:8;
        BN2002(i,:) = [nbpiid(X,rmax,i,0)   nbplog(X,rmax,i,0)];
    end;
    estim_numfac.bn2002 = BN2002;
    
    % Onatski (2010)
    Onatski2010 = onatski(X,rmax);
    estim_numfac.o2010 = Onatski2010;
    
    % Ahn-Horenstein (2013)
    [a,b]=Ahn_Horenstein(X,rmax);
    AH2013 = [a b];
    estim_numfac.ah2013 = AH2013;
    
    % Hallin-Liska (2007)
    HL2007 = HLestimate(X',rmax);
    estim_numfac.hl2007 = HL2007;
    
%     % Parallel test
%     ndatsets = 1000; percent = 95; info=0;
%     iter=0;
%     for ii=1:2
%         kind=ii;
%         for jj=1:2
%             iter=iter+1;
%             randtype = jj;
%             [numfac,realeval,means,percentiles]=parallel_test(X,ndatsets,percent,kind,randtype,info);
%             paral_nf(1,iter)=numfac;
%         end
%     end
%     % paral: (kind=1,rand=1),(kind=1,rand=2),(kind=2,rand=1),(kind=2,rand=2)
%     Parallel(:,tt) = paral_nf';
%     
%     % MAP test
%     [m1,m2] = map_test(corr(X));
%     map_nf(:,tt) = [m1,m2]';
%     % map_nf: (squared part corr), (^4 part corr)
    
    % %% Estimate number of factors ABC (2010)
    num_f = ABC_crit(X,rmax,1,floor(size(X,2)/4));
    estim_numfac.abc = num_f(end);
    
    % Bai-Ng (2007)
    [ehat,Fhat,lamhat,ve2]=pc(X,rmax); % get rhat factors   using L'L/N normalization
    [e,beta]=est_e(T,rmax-2,0,Fhat); % estimate VAR in Fht
    xepsilon=[2 2;1 1;.5 .5];     %covariance matrix
    eps=xepsilon/min([N^(.4);T^(.4)]);
    [stat1,stat2,evals_d]=dfactest(e,N,T,eps,1);
    BN2007 = [stat1' stat2'];
    estim_numfac.b2007 = BN2007;
    
    % Amengual-Watson (2007)
    [qhat]=dfactest_AW(X,rmax,rmax-2,1,1,0);
    AW2007= qhat;
    estim_numfac.aw2007 = AW2007;
    
    
    
    