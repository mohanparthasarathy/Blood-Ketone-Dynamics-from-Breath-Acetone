function study = run_inverse_study(outputDir,opts)
%RUN_INVERSE_STUDY Top-level entry point for synthetic inverse-study data.
if nargin<1 || isempty(outputDir)
 outputDir=fullfile(fileparts(mfilename('fullpath')),'outputs','cohort');
end
if nargin<2; opts=struct; end
cohort=generate_synthetic_cohort(outputDir,opts);
report=run_recoverability_analysis(cohort,outputDir);
study=struct('cohort',cohort,'recoverability',report);
end
