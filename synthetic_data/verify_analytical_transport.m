function report = verify_analytical_transport(makePlot)
%VERIFY_ANALYTICAL_TRANSPORT Compare MOL with the Ogata--Banks solution.
% PDE: c_s + c_x = Pe^{-1} c_xx on a semi-infinite domain, with c(0,s)=1,
% c(x,0)=0. We initialize at s0>0 from the exact solution to avoid testing a
% discontinuous corner. The truncated far boundary is also set to the exact value.
if nargin<1; makePlot=true; end
Pe=30; s0=0.02; sf=0.35; xMax=3; grids=[61 121 241 481];
err=zeros(size(grids)); dxs=err;
for g=1:numel(grids)
    N=grids(g); x=linspace(0,xMax,N).'; dx=x(2)-x(1); dxs(g)=dx;
    y0=exactOB(x(2:N-1),s0,Pe);
    opts=odeset('RelTol',1e-10,'AbsTol',1e-12);
    [~,Y]=ode15s(@rhs,[s0 sf],y0,opts);
    numeric=[1;Y(end,:).';exactOB(xMax,sf,Pe)];
    exact=exactOB(x,sf,Pe);
    interior=x<=1;
    err(g)=sqrt(mean((numeric(interior)-exact(interior)).^2));
end
orders=log(err(1:end-1)./err(2:end))./log(dxs(1:end-1)./dxs(2:end));
report=table(grids(:),dxs(:),err(:),[NaN;orders(:)], ...
    'VariableNames',{'N','dx','L2_error','observed_order'});
disp(report)
if makePlot
    figure('Color','w'); loglog(dxs,err,'o-','LineWidth',1.5); grid on
    xlabel('\Delta\xi'); ylabel('L_2 error'); title('Ogata--Banks verification');
end
    function dy=rhs(s,y)
        c=[1;y;exactOB(xMax,s,Pe)];
        dy=-(c(2:N-1)-c(1:N-2))/dx + ...
            (c(3:N)-2*c(2:N-1)+c(1:N-2))/(Pe*dx^2);
    end
end

function c=exactOB(x,s,Pe)
% Stable form of 0.5*(erfc(z1)+exp(Pe*x)*erfc(z2)).
z1=(x-s)./(2*sqrt(s/Pe)); z2=(x+s)./(2*sqrt(s/Pe));
c=0.5*(erfc(z1)+exp(-Pe*(x-s).^2./(4*s)).*erfcx(z2));
end
