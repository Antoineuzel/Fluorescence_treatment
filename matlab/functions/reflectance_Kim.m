function res=reflectance_Kim(x,rho,ua,lambda,n_medium,n_probe,wl_min,wl_max) 
    c=x(1:2);
    d=x(3);
    p=x(4:5);
    idx_min=find(min((lambda-wl_min).^2)==(lambda-wl_min).^2);
    idx_max=find(min((lambda-wl_max).^2)==(lambda-wl_max).^2);
    lambda=lambda(:,idx_min:idx_max);
    if length(n_probe(1)>1)
        n_probe=n_probe(:,idx_min:idx_max);
    end
    ua=ua(:,idx_min:idx_max);
    n=n_medium./n_probe;
    Rf=0.0636*n+0.668+0.710./n-1.44./(n.^2);
    kd=(1+Rf)./(1-Rf);
    uscat=c(1)*(d*(lambda(1,:)/600).^(-4)+(1-d)*(lambda(1,:)/600).^-c(2));
    ublood=p(1)*(p(2)*ua(1,:)+(1-p(2))*ua(2,:));% blood concentration and oxygen saturation
%     uabs=ublood+x(6)*ua(3,:)+x(7)*ua(4,:)+x(8)*ua(5,:);
    uabs=ublood+x(7)*ua(4,:)+x(8)*ua(5,:);
%     uabs=ublood;
    a=uscat./(uscat+uabs);
    z0=1./(uscat);
    D=1./(3*(uscat));
    zb=2*kd.*D;
    r1=sqrt(z0.^2+rho^2);
    r2=sqrt((z0+2*zb).^2+rho^2);
    ueff=sqrt(3*uabs.*uscat);
    kap=sqrt(3*uabs.*ueff);
    res1=z0.*(ueff+1./r1).*(exp(-ueff.*r1))./(r1.^2);
    res2=(z0+2*zb).*(ueff+1./r2).*(exp(-ueff.*r2))./(r2.^2);
    res=a/(4*pi).*(res1+res2);
end