function report = verify_manufactured_solution(makePlot)
%VERIFY_MANUFACTURED_SOLUTION MMS for full coupled gas--wall PDE/ODE.
% The forcing is constructed so c=1+0.2*sin(pi*x/2)*exp(-s) and
% w=0.8+0.1*cos(pi*x)*exp(-s/2) solve the forced equations exactly.
if nargin<1; makePlot=true; end
Pe=30; Daw=0.8; Dwa=0.8; Dwb=0.25; chi=0.1; a=1; sf=0.7;
grids=[31 61 121 241]; errC=zeros(size(grids)); errW=errC; dxs=errC;
for g=1:numel(grids)
    N=grids(g); x=linspace(0,1,N).'; dx=x(2)-x(1); dxs(g)=dx;
    [c0,w0]=exactFields(x,0); y0=[c0(2:N);w0];
    opts=odeset('RelTol',2e-10,'AbsTol',1e-12);
    [~,Y]=ode15s(@rhs,[0 sf],y0,opts);
    cn=[1;Y(end,1:N-1).']; wn=Y(end,N:2*N-1).';
    [ce,we]=exactFields(x,sf);
    errC(g)=sqrt(mean((cn-ce).^2)); errW(g)=sqrt(mean((wn-we).^2));
end
ordC=log(errC(1:end-1)./errC(2:end))./log(2);
ordW=log(errW(1:end-1)./errW(2:end))./log(2);
report=table(grids(:),dxs(:),errC(:),[NaN;ordC(:)],errW(:),[NaN;ordW(:)], ...
 'VariableNames',{'N','dx','L2_c','order_c','L2_w','order_w'});
disp(report)
if makePlot
    figure('Color','w'); loglog(dxs,errC,'o-',dxs,errW,'s-','LineWidth',1.5);
    grid on; xlabel('\Delta\xi'); ylabel('L_2 error');
    title('Manufactured-solution convergence'); legend('gas c','wall w');
end
    function dy=rhs(s,y)
        c=[1;y(1:N-1)]; w=y(N:2*N-1);
        dc=zeros(N,1);
        dc(2:N-1)=-(c(2:N-1)-c(1:N-2))/dx + ...
          (c(3:N)-2*c(2:N-1)+c(1:N-2))/(Pe*dx^2)-Daw*(c(2:N-1)-w(2:N-1));
        dc(N)=-(c(N)-c(N-1))/dx+2*(c(N-1)-c(N))/(Pe*dx^2)-Daw*(c(N)-w(N));
        dw=Dwa*(c-w)-Dwb*(w-chi*a);
        [fc,fw]=forcing(x,s);
        dy=[dc(2:N)+fc(2:N);dw+fw];
    end
    function [c,w]=exactFields(xv,s)
        c=1+0.2*sin(pi*xv/2)*exp(-s);
        w=0.8+0.1*cos(pi*xv)*exp(-s/2);
    end
    function [fc,fw]=forcing(xv,s)
        [c,w]=exactFields(xv,s);
        cs=-0.2*sin(pi*xv/2)*exp(-s);
        cx=0.1*pi*cos(pi*xv/2)*exp(-s);
        cxx=-0.05*pi^2*sin(pi*xv/2)*exp(-s);
        ws=-0.05*cos(pi*xv)*exp(-s/2);
        fc=cs+cx-cxx/Pe+Daw*(c-w);
        fw=ws-Dwa*(c-w)+Dwb*(w-chi*a);
    end
end
