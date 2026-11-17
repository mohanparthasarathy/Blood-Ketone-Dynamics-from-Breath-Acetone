function out = solve_fast_breath(a_n,p,U)
%SOLVE_FAST_BREATH Method-of-lines airway solve for one frozen slow snapshot.
% U is dimensionless exhalation velocity u/u_star. It may be a scalar or a
% function handle U(s). The default is U=1.
if nargin < 3 || isempty(U); U = 1; end
N = p.fast.nx;
xi = linspace(0,1,N).'; dx = xi(2)-xi(1);
sEnd = p.fast.breath_duration_s/p.scales.T_p_s;
sEval = linspace(0,sEnd,301).';

% y = [c(2:N); w(1:N)], while c(1)=a_n is imposed exactly.
y0 = zeros(2*N-1,1);
opts = odeset('RelTol',2e-7,'AbsTol',1e-9,'NonNegative',1:(2*N-1));
[s,y] = ode15s(@rhs,sEval,y0,opts);
c = [a_n*ones(numel(s),1), y(:,1:N-1)].';
w = y(:,N:2*N-1).';

out.xi = xi; out.x_m = xi*p.scales.L_m;
out.s = s; out.time_s = s*p.scales.T_p_s;
out.c = c; out.w = w;
out.C_ppm = c*p.scales.C_star_ppm;
out.C_exh_ppm = out.C_ppm(end,:).';

    function dy = rhs(sNow,yv)
        if isa(U,'function_handle'); uNow=U(sNow); else; uNow=U; end
        if ~isscalar(uNow) || ~isfinite(uNow) || uNow<=0
            error('U(s) must be a finite positive scalar.');
        end
        cNow = [a_n; yv(1:N-1)];
        wNow = yv(N:2*N-1);
        dc = zeros(N,1);
        dc(2:N-1) = -uNow*(cNow(2:N-1)-cNow(1:N-2))/dx + ...
            (cNow(3:N)-2*cNow(2:N-1)+cNow(1:N-2))/(p.fast.Pe*dx^2) - ...
            p.fast.Da_aw*(cNow(2:N-1)-wNow(2:N-1));
        % Reflected ghost cell enforces dc/dxi=0 at xi=1.
        dc(N) = -uNow*(cNow(N)-cNow(N-1))/dx + ...
            2*(cNow(N-1)-cNow(N))/(p.fast.Pe*dx^2) - ...
            p.fast.Da_aw*(cNow(N)-wNow(N));
        dw = p.fast.Da_wa*(cNow-wNow)- ...
             p.fast.Da_wb*(wNow-p.fast.chi*a_n);
        dy = [dc(2:N); dw];
    end
end
