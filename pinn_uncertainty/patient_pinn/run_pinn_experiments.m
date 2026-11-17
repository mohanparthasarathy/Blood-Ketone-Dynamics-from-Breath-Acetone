%RUN_PINN_EXPERIMENTS Oracle validation followed by practical 10% prior.
clear; clc
root=fileparts(mfilename('fullpath'));

oracleOpts=struct('gamma_mode',"oracle",'max_iterations',8000, ...
    'num_restarts',3,'refit_iterations',1000,'use_lbfgs',true, ...
    'num_collocation',161,'fast_nx',121,'print_every',500);
oracle_pinn=run_patient_pinn([],fullfile(root,'outputs','oracle'),oracleOpts);

priorOpts=struct('gamma_mode',"prior",'gamma_prior_cv',.10, ...
    'max_iterations',8000,'num_restarts',3,'refit_iterations',1000, ...
    'use_lbfgs',true,'num_collocation',161,'fast_nx',121,'print_every',500);
prior10_pinn=run_patient_pinn([],fullfile(root,'outputs','prior10'),priorOpts);

comparison=[oracle_pinn.summary;prior10_pinn.summary];
writetable(comparison,fullfile(root,'outputs','pinn_comparison.csv'));
disp(comparison)
save(fullfile(root,'outputs','pinn_experiments.mat'), ...
    'oracle_pinn','prior10_pinn','comparison','-v7.3');
