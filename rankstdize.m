function rankvar = rankstdize(rawvar)

[T,N]       = size(rawvar);
rankvar     = nan(T,N);
for t=1:T
    n               = sum(~isnan(rawvar(t,:)));

    rk              = tiedrank(rawvar(t,:));
    rankvar(t,:)    = rk/n - nanmean(rk/n);

end
    

    