function final_smoke=run_final_smoke_test
%RUN_FINAL_SMOKE_TEST Fast structural test before the full resumable run.
here=fileparts(mfilename('fullpath'));
o=struct('mode',"smoke", ...
    'output_dir',fullfile(here,'outputs','smoke'), ...
    'resume',true);
final_smoke=run_remaining_project(o);
end

