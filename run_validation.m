function run_validation()
%RUN_VALIDATION Verify the forward model and compare direction rules.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'synthetic_data'),'-begin');
run_verification_suite;
cohort=fullfile(root,'synthetic_data','outputs','cohort','synthetic_cohort.mat');
run_recoverability_analysis(cohort,fullfile(root,'synthetic_data','outputs','cohort'));
end
