function [state,dstate]=evaluate_pinn(net,time_h,tmax)
%EVALUATE_PINN Return states and their time derivatives.
tau=dlarray(single(time_h(:).'),'CB');
[Y,dY]=dlfeval(@differentiate,net,tau,tmax);
state=double(gather(extractdata(Y)));
dstate=double(gather(extractdata(dY)));
end
function [Y,dY]=differentiate(net,tau,tmax)
X=stripdims(tau); Y=stripdims(forward(net,2*tau/tmax-1));
d1=reshape(dljacobian(Y(1,:),X,1,EnableHigherDerivatives=false),1,[]);
d2=reshape(dljacobian(Y(2,:),X,1,EnableHigherDerivatives=false),1,[]);
d3=reshape(dljacobian(Y(3,:),X,1,EnableHigherDerivatives=false),1,[]);
dY=[d1;d2;d3];
end
