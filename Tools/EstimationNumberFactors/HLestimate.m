%cd 'Francisco'
% This program loads the data simulated from Bouakez, Cardia and
% Ruge-Murcia's model and applies the Hallin-Liska
% procedure to determine the number of dynamic factors in these data
% The procedure is as described in Hallin and Liska (2007) "Determining
% the number of factors in the general dynamic factor model", Journal
% of the American Statistical Association, 603-617 (HL henceforth)


function [out] = HLestimate(dataclean,qmax);

%% Output:
% out - estimated number of factors

%% Inputs:
% dataclean - N x T matrix of data !!! ATTENTION, LES DIMENSIONS SONT N X T !!!
% qmax - nombre maximal de facteurs


[nn,TT]=size(dataclean);


c=0.01:0.01:3;
[cn,cm]=size(c);

% Consider subsamples of the data (as described in point 4 on page 610
% of HL (2007))
for jj=3:-1:0
    n=nn-jj*10;
    T=TT-jj*10;
    data=dataclean(1:n,1:T);
    % Choose MT as described in point 2 on page 610)
    MT=floor(0.7*sqrt(T));
    % Choose the grid c, and a penalty as described in point 3
    % on page 610
    p1=(MT^(-2)+sqrt(MT)/sqrt(T)+1/n)*log(min([n;MT^2;sqrt(T)/sqrt(MT)]));
    % compute the sum of squared residuals from fitting different
    % number of dynamic factors, average over frequency grid
    % consider frequency grid
    GAMMA=[];
    for u=0:MT
        GAMMA=[GAMMA;data(:,u+1:end)*data(:,1:end-u)'/(2*pi*T)];
    end
    LAMBDA=[];
    for l=0:MT
        theta=pi*l/(MT+0.5);
        vect=(1-(1:MT)/MT).*exp(-sqrt(-1)*(1:MT)*theta);
        vect=kron(vect,eye(n));
        SIGMA1=vect*GAMMA(n+1:end,:);
        SIGMA=GAMMA(1:n,:)+SIGMA1+SIGMA1';
        % compute eigenvalues of the spectral density matrix
        lambda=eig(SIGMA);
        % make sure they are real (computer may do some numerical error and give you
        % a number with small imaginary part)
        lambda=real(lambda);
        % sort the eigenvalues
        lambdas=sort(lambda);
        % reorder so that the largest eigenvalues come first
        LAMBDA(:,l+1)=flipdim(lambdas,1);
    end
    for k=0:qmax
        % compute sum of squared residuals (do the cicle to average
        % over the frequency grid)
        IC1(k+1,1)=(sum(LAMBDA(k+1:end,1))+2*sum(sum(LAMBDA(k+1:end,2:end))))/(n*(2*MT+1));
    end
    % compute the information criteria corresponding to
    % different scales on the grid 1:cm
    for ii=1:cm
        IC21(:,ii)=log(IC1(:,1))+(0:qmax)'*c(1,ii)*p1;
    end
    % find the number of factors (minimum of the criteria) corresponding
    % to different scales
    [waste,q21]=min(IC21);
    Q21(jj+1,:)=q21-1;
end
%Solution to ASSIGNMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do averaging described in equation (10) on page 607
S221=sum((Q21-ones(4,1)*sum(Q21)/4).^2)/4;
%END OF solution to ASSIGNMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now, do the second stability region procedure
iii=1;
ff=S221(1,iii);
while ff==0
    if iii==cm
        ff=1;
    else
        iii=iii+1;
        ff=S221(1,iii);
    end
end
while ff~=0
    if iii==cm
        ff=0;
    else
        iii=iii+1;
        ff=S221(1,iii);
    end
end
out=Q21(1,iii);
%disp('number of dynamic factors is');
%disp(out)

