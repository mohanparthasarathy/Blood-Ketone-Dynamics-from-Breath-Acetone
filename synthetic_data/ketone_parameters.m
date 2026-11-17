function p = ketone_parameters()
%KETONE_PARAMETERS Baseline scales and dimensionless parameters.
p.scales.K_star_mM = 1.0;
p.scales.T_m_h = 1.0;
p.scales.L_m = 0.30;
p.scales.u_star_m_s = 0.50;
p.scales.lambda_ba = 250.0;
p.scales.optical_path_m = 0.10;

p.slow.rho_B = 0.80;
p.slow.kappa_BQ = 0.15;
p.slow.kappa_QB = 0.65;
p.slow.kappa_Q = 0.35;
p.slow.kappa_QA = 0.30;
p.slow.Gamma_A = 1.50;
p.slow.r0 = 0.22;
p.slow.delta_r = 2.10;
p.slow.tau_c = 3.0;
p.slow.k_r = 2.2;

% A_star is derived, not independently selected.
p.scales.A_star_mM = p.slow.kappa_QA*p.scales.K_star_mM/p.slow.Gamma_A;
p.scales.T_p_s = p.scales.L_m/p.scales.u_star_m_s;
p.scales.C_star_ppm = p.scales.A_star_mM/p.scales.lambda_ba*24450;
p.scales.epsilon_t = p.scales.T_p_s/(3600*p.scales.T_m_h);

p.fast.Pe = 30.0;
p.fast.Da_aw = 0.80;
p.fast.Da_wa = 0.80;
p.fast.Da_wb = 0.25;
p.fast.chi = 0.10;
p.fast.nx = 121;
p.fast.breath_duration_s = 5.0;
p.fast.D_m2_s = p.scales.u_star_m_s*p.scales.L_m/p.fast.Pe;
end
