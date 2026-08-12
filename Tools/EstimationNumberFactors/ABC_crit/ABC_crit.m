% Modification of the criterion by Bai and Ng (2002) for determining number of static factors 
% AUTHORS: Lucia Alessi, Matteo Barigozzi, Marco Capasso
% LAST UPDATE: June, 25 2009

function [num_f]=ABC_crit(panel,kmax,cmax,nbck)

% INPUT     panel: (T x n) stationary panel 
%           kmax : upper bound on the number of factors, e.g. kmax = 10  
%           cmax : c = [0:cmax], e.g cmax = 5
%           nbck : T x n_j subpanels are used where
%                  n_j = n - nbck : n
% OUTPUT    num_f(end,:) is the number of factors chosen by the IC1 criterion 
%           for different values of c (red solid line in the figures) and
%           n_j = N
%           std(num_f) is the variance in the number of factors chosen by the
%           IC1 criterion for different values of c when varying the size of the
%           subpanels i.e. for n_j -> N (blue dashed line)

[T,n] = size(panel);
x = (panel - ones(T,1)*mean(panel))./(ones(T,1)*std(panel,1));

s=0;

for Nsub = n-nbck:n
    disp(sprintf('subsample size %d',Nsub));
    s = s+1;
    [a Ns]=sort(rand(n,1));
    subpanel = x(1:T,Ns(1:Nsub));
    m_subpanel = mean(subpanel);
    s_subpanel = std(subpanel);
    subpanel = (subpanel - ones(T,1)*m_subpanel)./(ones(T,1)*s_subpanel);   % standardize the subpanel

    V = subr_1(subpanel,kmax); 
    [a1,a2] = size(V); 
  
    penalty = ((Nsub+T)/(Nsub*T))*log((Nsub*T)/(Nsub+T));
    
    Vmax = V(kmax,1);
    
    for c = 1:floor(cmax*100)
        cc = c/100;
          
        IC_1 =zeros(kmax+1,1); % when the tested number of factors is zero, log(V)=0 and k*p*c=0
        IC_1(2:kmax+1,1) = log(V) + kron([1:kmax]',ones(1,a2)).*penalty*cc;   % for IC the penalty is k*p*c
        
        [rr,rs]=find((IC_1 == ones(kmax+1,1)*min(IC_1))==1);
        num_f(s,c)=rr-1;    
    end 
end 

set(0,'DefaultLineLineWidth',2);
cr=[1:floor(cmax*100)]'/100;

% figure 
% plot(cr,num_f(end,:),'r-')
% axis tight
% hold all
% plot(cr,5*std(num_f),'b--')
% xlabel('c')
% axis tight
% legend('r^{*T}_{c;N}','S_c')
% title('ABC estimated number of factors')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  V = subr_1(x,kmax)
[T,Nsub] = size(x);

for k=1:kmax
    GG=cov(x);
    [autovet,autoval]=eig(GG);
    [autoval,IXX] = sort((diag(autoval)));   % sort eigenvalues and eigenvectors
    autoval = flipud(autoval);
    IXX = flipud(IXX);
    autovet = autovet(:,IXX);                % now arrange the eigenvalues  and eigenvectors according to the decreasing order of eigenvalues
    autovet = autovet(:,1:k);
    lambda = autovet*sqrt(Nsub);
    clear autoval autovet IXX
    F = (x*lambda)/Nsub;                          % see lines 7-12 of page 9 (198) of Bai-Ng 2002
    newF = F*((F'*F)/T)^(1/2);
    newlambda = (newF'*newF)\newF'*x;             % see line 4 of page 10 (199) of Bai_Ng 2002
    newlambda=newlambda';
    clear F lambda
    SQ=0;
    for i=1:Nsub
        for t=1:T
            SQ=SQ+(x(t,i)-newlambda(i,:)*(newF(t,:))')^2;
        end
    end
    clear newlambda f
    V(k,1)=(SQ)/(Nsub*T);   
end