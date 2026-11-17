function out = solve_slow_metabolism(p, scenario)
%SOLVE_SLOW_METABOLISM Solve dimensionless BHB/AcAc/acetone dynamics.
if nargin < 2; scenario = 'developing'; end
tau = linspace(0,8,481).';
y0 = [0.18; 0.06; 0.06];
opts = odeset('RelTol',1e-9,'AbsTol',1e-11,'NonNegative',1:3);
[tau,y] = ode15s(@rhs,tau,y0,opts);

out.tau = tau;
out.time_h = tau*p.scales.T_m_h;
out.r = ketogenic_drive(tau,p.slow,scenario);
out.b = y(:,1); out.q = y(:,2); out.a = y(:,3);
out.B_mM = p.scales.K_star_mM*out.b;
out.Q_mM = p.scales.K_star_mM*out.q;
out.A_blood_mM = p.scales.A_star_mM*out.a;

    function dy = rhs(t,yv)
        r = ketogenic_drive(t,p.slow,scenario);
        b = yv(1); q = yv(2); a = yv(3);
        dy = [p.slow.rho_B*r-(1+p.slow.kappa_BQ)*b+p.slow.kappa_QB*q; ...
              (1-p.slow.rho_B)*r+p.slow.kappa_BQ*b- ...
              (p.slow.kappa_QB+p.slow.kappa_Q+p.slow.kappa_QA)*q; ...
              p.slow.Gamma_A*(q-a)];
    end
end

function r = ketogenic_drive(tau,p,scenario)
z = max(-60,min(60,p.k_r*(tau-p.tau_c)));
switch lower(scenario)
    case 'developing'
        r = p.r0+p.delta_r./(1+exp(-z));
    case 'resolving'
        r = p.r0+p.delta_r./(1+exp(z));
    case 'baseline'
        r = p.r0+zeros(size(tau));
    otherwise
        error('Unknown scenario: %s',scenario)
end
end
