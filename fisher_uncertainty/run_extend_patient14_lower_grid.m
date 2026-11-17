% Extend only the lower side of the refined patient-14 profile.
% The 11 completed points are discovered and loaded from their checkpoints.
clear; clc
here=fileparts(mfilename('fullpath')); projectDir=fileparts(here);
R=load(fullfile(projectDir,'fisher_pinn_validation','outputs', ...
    'fisher_validation_full','fisher_pinn_validation.mat'));
j=find(string(R.validation.summary.design)=="fisher_D8",1);
F=R.validation.results{j};
q=F.parameters(F.parameters.patient_id==14 & F.parameters.parameter=="Gamma_A",:);
center=q.estimate(1);
completedGrid=center*linspace(0.45,1.55,11);
extendedGrid=sort(unique([0.30 0.40 0.50 0.60 completedGrid]));

opts=struct('patient_ids',14,'gamma_grid_values',extendedGrid, ...
    'profile_iterations',5000,'profile_min_iterations',2500, ...
    'profile_restarts',2,'profile_refit_iterations',500, ...
    'profile_use_lbfgs',true,'profile_lbfgs_iterations',150, ...
    'num_collocation',121,'fast_nx',81,'print_every',1000,'resume',true);

fisher_uncertainty_extended=run_fisher_uncertainty(F,[], ...
    fullfile(here,'outputs','refined_patient14'),opts);
disp(fisher_uncertainty_extended.summary)
disp(fisher_uncertainty_extended.coverage)
disp(fisher_uncertainty_extended.profile)
disp(fisher_uncertainty_extended.profile( ...
    fisher_uncertainty_extended.profile.accepted,:))
disp(fisher_uncertainty_extended.alerts)
