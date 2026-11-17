function report=check_patient_pinn_setup()
%CHECK_PATIENT_PINN_SETUP Verify dljacobian layout against finite differences.
assert(exist('dljacobian','file')==2, ...
    'Requires MATLAB R2024b+ and Deep Learning Toolbox.');
net=dlnetwork([featureInputLayer(1,Normalization="none"); ...
    fullyConnectedLayer(8);tanhLayer;fullyConnectedLayer(3)]);
t=linspace(.2,7.8,9); [~,ad]=evaluate_pinn(net,t,8);
h=1e-3; yp=extractdata(forward(net,dlarray(single(2*(t+h)/8-1),'CB')));
ym=extractdata(forward(net,dlarray(single(2*(t-h)/8-1),'CB')));
fd=double((yp-ym)/(2*h)); err=max(abs(ad-fd),[],'all');
report=table(err,err<5e-3,'VariableNames',{'max_derivative_error','passed'});
disp(report)
assert(report.passed,'Automatic-derivative smoke test failed.');
end
