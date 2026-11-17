function out=decode_pinn_spectra(noisy,clean,wn,p)
%DECODE_PINN_SPECTRA Known one-VOC spectral inversion; clean is scoring only.
wn=double(wn(:).');
shape=.55*exp(-.5*((wn-1215)/30).^2)+.38*exp(-.5*((wn-1365)/38).^2)+ ...
    exp(-.5*((wn-1730)/24).^2);
g=linspace(1000,1900,451);
full=.55*exp(-.5*((g-1215)/30).^2)+.38*exp(-.5*((g-1365)/38).^2)+ ...
    exp(-.5*((g-1730)/24).^2);
theta=8e-5*(p.scales.optical_path_m/.10)*shape/max(full); den=sum(theta.^2);
nT=size(noisy,1); n=size(noisy,3); Chat=zeros(nT,n); Ctrue=Chat;
for j=1:n
    Chat(:,j)=max(0,double(noisy(:,:,j))*theta.'/den);
    Ctrue(:,j)=max(0,double(clean(:,:,j))*theta.'/den);
end
out=struct('C_hat_ppm',Chat,'C_true_ppm',Ctrue);
end
