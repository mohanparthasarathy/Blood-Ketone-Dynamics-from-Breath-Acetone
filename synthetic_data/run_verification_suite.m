function verification = run_verification_suite(outputDir)
%RUN_VERIFICATION_SUITE Run analytical and manufactured-solution tests.
if nargin<1
    outputDir=fullfile(fileparts(mfilename('fullpath')),'outputs');
end
if ~exist(outputDir,'dir'); mkdir(outputDir); end
verification.slow_exact=verify_analytical_slow();
verification.analytical=verify_analytical_transport(false);
verification.mms=verify_manufactured_solution(false);
writetable(verification.slow_exact,fullfile(outputDir,'slow_exact_errors.csv'));
writetable(verification.analytical,fullfile(outputDir,'analytical_convergence.csv'));
writetable(verification.mms,fullfile(outputDir,'mms_convergence.csv'));
save(fullfile(outputDir,'verification_results.mat'),'verification');
end
