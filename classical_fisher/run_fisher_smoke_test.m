%RUN_FISHER_SMOKE_TEST Representative patient-14 classical design.
clear; clc
here=fileparts(mfilename('fullpath')); root=fileparts(here);
opts=struct('patient_id',14,'fast_nx',41,'phase_candidates',11, ...
    'wavelength_candidates',31,'wavelength_subset_sizes',[3 5 8], ...
    'time_subset_sizes',[2 3 5],'phase_subset_sizes',[3 5 8]);
fisher_smoke=run_classical_fisher_design([], ...
    fullfile(here,'outputs','smoke_patient14'),opts);
disp(fisher_smoke.diagnostics)
disp(fisher_smoke.wavelength_designs)
disp(fisher_smoke.time_design)
disp(fisher_smoke.phase_design)
