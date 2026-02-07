function [ret_total,ret_mean,SR,all_L,all_weight,Ftil_total,Ftil_mean,SR_Ftil,IR,alpha,pval,stats]=main(retd,fffac,rf)

rollwin = 120;



[fullT,N] = size(retd);

momwin      = 20;
fwdwin      = 20;


Tno     = floor(fullT/fwdwin);
fwdix   = kron((1:Tno)',ones(fwdwin,1));
fwdix   = [fwdix;nan(fullT-length(fwdix),1)];
Rfwd    = nan(Tno,N);   
S       = nan(Tno,N);

FFfwd   = nan(Tno,size(fffac,2));
skip    = 2;


for t=3:Tno
    loc         = find(fwdix==t);
    Rfwd(t,:)   = prod(1+retd(loc,:),1)-1;
    sigloc      = (loc(1)-skip-(momwin-1)):(loc(1)-skip);
    S(t,:)      = sum(retd(sigloc,:),1);
    
    mkt         = fffac(loc,1)+rf(loc);
    ffoth       = fffac(loc,2:end);
    FFfwd(t,:)  = prod(1+[mkt ffoth],1)-1;
    RFfwd       = prod(1+rf(loc),1)-1;
    FFfwd(t,1)  = FFfwd(t,1)-RFfwd;
    
end


all_weight=zeros(Tno-rollwin-2,N);
all_L=zeros(Tno-rollwin-2,N,N);
ret_total=zeros(Tno-rollwin-2,1);
Ftil_total=zeros(Tno-rollwin-2,1);
FFfwd_total=FFfwd(rollwin+3:Tno,:);


for tau = rollwin+2:Tno-1
    Rtrnslc = Rfwd(tau-rollwin+1:tau,:); 
    Strnslc = S(tau-rollwin+1:tau,:);
    Rtstslc = Rfwd(tau+1,:);
    Ststslc = S(tau+1,:);
    
    Strnslc            = rankstdize(Strnslc);
    Ststslc            = rankstdize(Ststslc);
    
    Rtrnslc            = Rtrnslc - nanmean(Rtrnslc,2); 
    Rtstslc            = Rtstslc - nanmean(Rtstslc,2); 
    
    Pie_hat=Strnslc'*Rtrnslc/rollwin;
    [L,~]=fixedpoint(Pie_hat);
    all_L(tau-rollwin-1,:,:)=L;
    weight=Ststslc*L';
    all_weight(tau-rollwin-1,:)=weight;
    ret_new=weight*Rtstslc';
    ret_total(tau-rollwin-1)=ret_new;
    Ftil_total(tau-rollwin-1)=Ststslc*Rtstslc';
end



ret_mean=mean(ret_total);
SR=ret_mean/std(ret_total);

Ftil_mean=mean(Ftil_total);
SR_Ftil=Ftil_mean/std(Ftil_total);

stats= regstats(ret_total,[Ftil_total,FFfwd_total],'linear',{'rsquare','tstat','r','fstat'});
IR= sharpe(stats.tstat.beta(1)+stats.r);
alpha= stats.tstat.beta(1);
pval=stats.tstat.pval(1)/2;




