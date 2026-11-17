function run_reconstruction()
%RUN_RECONSTRUCTION Fit the two response-calibration conditions.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'synthetic_data'));
addpath(fullfile(root,'patient_pinn'),'-begin');
cohort=fullfile(root,'synthetic_data','outputs','cohort','synthetic_cohort.mat');
settings=struct('max_iterations',8000,'num_restarts',3, ...
    'refit_iterations',1000,'use_lbfgs',true);
settings.gamma_mode="oracle";
run_patient_pinn(cohort,fullfile(root,'patient_pinn','outputs','oracle_v2'),settings);
settings.gamma_mode="prior"; settings.gamma_prior_cv=0.10;
run_patient_pinn(cohort,fullfile(root,'patient_pinn','outputs','prior10_v2'),settings);
end
