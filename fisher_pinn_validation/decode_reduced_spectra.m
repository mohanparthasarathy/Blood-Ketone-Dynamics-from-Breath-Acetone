function out=decode_reduced_spectra(noisy,clean,wn,p,fitNuisance)
%DECODE_REDUCED_SPECTRA Spectral inversion with optional nuisance removal.
if nargin<5; fitNuisance=false; end
wn=double(wn(:).');
shape=.55*exp(-.5*((wn-1215)/30).^2)+.38*exp(-.5*((wn-1365)/38).^2)+ ...
    exp(-.5*((wn-1730)/24).^2);
g=linspace(1000,1900,451);
full=.55*exp(-.5*((g-1215)/30).^2)+.38*exp(-.5*((g-1365)/38).^2)+ ...
    exp(-.5*((g-1730)/24).^2);
theta=8e-5*(p.scales.optical_path_m/.10)*shape/max(full);
if fitNuisance
    x=(wn-mean(wn))/max(eps,max(wn)-min(wn));
    H=[ones(numel(wn),1),x(:),exp(-.5*((wn(:)-1640)/95).^2), ...
        exp(-.5*((wn(:)-1460)/42).^2)];
    H=H./max(abs(H),[],1); X=[theta(:),H]; Xpinv=pinv(X);
else
    den=sum(theta.^2);
end
nT=size(noisy,1); n=size(noisy,3); Chat=zeros(nT,n); Ctrue=Chat;
for j=1:n
    if fitNuisance
        beta=(Xpinv*double(noisy(:,:,j)).').';
        betaClean=(Xpinv*double(clean(:,:,j)).').';
        Chat(:,j)=max(0,beta(:,1)); Ctrue(:,j)=max(0,betaClean(:,1));
    else
        Chat(:,j)=max(0,double(noisy(:,:,j))*theta.'/den);
        Ctrue(:,j)=max(0,double(clean(:,:,j))*theta.'/den);
    end
end
out=struct('C_hat_ppm',Chat,'C_true_ppm',Ctrue);
end
