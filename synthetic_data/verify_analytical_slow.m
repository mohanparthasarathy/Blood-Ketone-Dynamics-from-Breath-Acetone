function report = verify_analytical_slow()
%VERIFY_ANALYTICAL_SLOW Compare ode15s against a matrix-exponential solution.
% For constant ketogenesis r0, y'=M*y+d has the exact solution
% y(t)=y_ss+expm(M*t)*(y0-y_ss).
p=ketone_parameters(); q=p.slow; r=q.r0;
M=[-(1+q.kappa_BQ), q.kappa_QB, 0; ...
    q.kappa_BQ, -(q.kappa_QB+q.kappa_Q+q.kappa_QA), 0; ...
    0, q.Gamma_A, -q.Gamma_A];
d=[q.rho_B*r; (1-q.rho_B)*r; 0];
y0=[0.18;0.06;0.06]; yss=-M\d;
t=linspace(0,8,161).'; exact=zeros(numel(t),3);
for k=1:numel(t)
    exact(k,:)=(yss+expm(M*t(k))*(y0-yss)).';
end
opts=odeset('RelTol',1e-10,'AbsTol',1e-12);
[~,numeric]=ode15s(@(~,y) M*y+d,t,y0,opts);
maxError=max(abs(numeric-exact),[],1);
report=table({'BHB';'AcAc';'acetone'},maxError.', ...
    'VariableNames',{'state','max_abs_error'});
disp(report)
end
