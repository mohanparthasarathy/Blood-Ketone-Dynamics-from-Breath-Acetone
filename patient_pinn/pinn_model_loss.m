function [loss,gradNet,gradRaw,parts,physical]=pinn_model_loss( ...
    net,raw,tauC,tauO,Cobs,gain,base,o,lb,ub,trueTheta,gammaPrior,scenario,w)
%PINN_MODEL_LOSS Curriculum-weighted data and governing-equation loss.
p=lb+(ub-lb).*sigmoid(raw); mode=string(o.gamma_mode);
if mode=="oracle"; gamma=dlarray(single(trueTheta(8)));
elseif mode=="population"; gamma=dlarray(single(o.gamma_population));
else; gamma=p(5); end
rpre=p(1); rpost=p(2); tc=p(3); kr=p(4);
if string(scenario)=="baseline"; rpost=rpre; end
physical=[rpre rpost tc kr gamma];

X=stripdims(tauC); Y=stripdims(forward(net,2*tauC/o.t_max_h-1));
b=Y(1,:); q=Y(2,:); a=Y(3,:);
db=reshape(dljacobian(b,X,1),1,[]); dq=reshape(dljacobian(q,X,1),1,[]);
da=reshape(dljacobian(a,X,1),1,[]);
z=max(-60,min(60,kr*(X-tc))); r=rpre+(rpost-rpre)./(1+exp(-z));
fB=db-(base.slow.rho_B*r-(1+base.slow.kappa_BQ)*b+base.slow.kappa_QB*q);
fQ=dq-((1-base.slow.rho_B)*r+base.slow.kappa_BQ*b- ...
    (base.slow.kappa_QB+base.slow.kappa_Q+base.slow.kappa_QA)*q);
fA=da-gamma*(q-a); physicsLoss=mean(fB.^2+fQ.^2+fA.^2,'all');

YO=forward(net,2*tauO/o.t_max_h-1); aO=YO(3,:);
Cstar=base.slow.kappa_QA*base.scales.K_star_mM/gamma/base.scales.lambda_ba*24450;
Cpred=gain.*Cstar.*aO; scale=max(.25,max(Cobs,[],'all')-min(Cobs,[],'all'));
dataLoss=mean(((Cpred-Cobs)/scale).^2,'all');
Y0=forward(net,dlarray(single(-1),'CB'));
icLoss=mean((Y0-dlarray(single(o.ic_target(:)))).^2,'all');
positiveLoss=mean(max(-Y,0).^2,'all'); gammaLoss=dlarray(single(0));
if mode=="prior"
    sigmaLog=sqrt(log(1+o.gamma_prior_cv^2));
    gammaLoss=((log(gamma)-log(gammaPrior))/sigmaLog).^2;
end
if string(scenario)=="developing"
    forcingLoss=((rpre-o.r_low_prior_center)/o.r_prior_scale).^2+ ...
        ((rpost-o.r_high_prior_center)/o.r_prior_scale).^2;
elseif string(scenario)=="resolving"
    forcingLoss=((rpre-o.r_high_prior_center)/o.r_prior_scale).^2+ ...
        ((rpost-o.r_low_prior_center)/o.r_prior_scale).^2;
else
    forcingLoss=((rpre-o.r_low_prior_center)/o.r_prior_scale).^2;
end
if string(scenario)~="baseline"
    forcingLoss=forcingLoss+((kr-o.kr_prior_center)/o.kr_prior_scale).^2;
end
loss=w.data*dataLoss+w.physics*physicsLoss+w.ic*icLoss+ ...
    w.positive*positiveLoss+w.gamma*gammaLoss+w.forcing*forcingLoss;
parts=[dataLoss physicsLoss icLoss gammaLoss];
[gradNet,gradRaw]=dlgradient(loss,net.Learnables,raw);
end
function y=sigmoid(x); y=1./(1+exp(-x)); end
