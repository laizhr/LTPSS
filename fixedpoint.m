function [L,Ltotal]=fixedpoint(Pie_hat)

kk=3;
theta=0.9999;



N=size(Pie_hat,1);

beta=100;
%gamma=0.1;

%gamma=0.1/25*N;

[U,lbd_pie,V] = svd(Pie_hat);
L=V(:,1:kk)*U(:,1:kk)';

eta=0.001;
gamma=eta*beta;


MaxIter=1000;
tol=1e-5;



k1 = 1;
RE1 = zeros(MaxIter+1,1);
RE1(1) = inf;
Ltotal = [];
while k1<=MaxIter && RE1(k1)>tol
            L_pre = L;

            L1=L+beta*Pie_hat';
            [U,lbd,V] = svd(L1);
            lbd=wthresh(diag(lbd),"s",gamma);
          %  L2=U*lbd*V';
            comp=[lbd,ones(N,1)];
            lbd=min(comp,[],2);
            
        
            L_fr=U*diag(lbd)*V';
            
            L=theta*L_fr+(1-theta)*L_pre;
     
            RE1(k1+1) = norm(L-L_pre,2)/norm(L_pre,2);
            Ltotal(k1,:,:) = L;
            k1 = k1+1;
end



end


