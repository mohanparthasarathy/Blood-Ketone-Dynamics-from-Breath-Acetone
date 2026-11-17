function validation=run_fisher_pinn_validation(outputDir,userOpts)
%RUN_FISHER_PINN_VALIDATION Compare full and reduced classical designs.
% All reduced designs use identical breath times and within-breath phases.

here=fileparts(mfilename('fullpath'));
addpath(here,'-begin');
if nargin<1 || isempty(outputDir)
    outputDir=fullfile(here,'outputs','fisher_validation');
end
if nargin<2; userOpts=struct; end
if ~exist(outputDir,'dir'); mkdir(outputDir); end

times=[0.5 1.25 2 2.75 3.5 5 5.75 8];
phases=[1.35 1.85 2.50 3.15 3.35 3.50 3.65 3.85 ...
    4.00 4.15 4.35 4.50 4.65 4.85 5.00];

designs=struct( ...
    'label',{"full_grid","fisher_D8","uniform_8","peak_8"}, ...
    'wn',{[], ...
        [1210 1226 1450 1630 1720 1736 1750 1900], ...
        [1000 1136 1256 1390 1510 1646 1766 1900], ...
        [1736 1720 1750 1706 1210 1226 1196 1240]}, ...
    'times',{[],times,times,times}, ...
    'phases',{[],phases,phases,phases});

allResults=cell(numel(designs),1); comparison=[]; byPatient=[];
for k=1:numel(designs)
    label=string(designs(k).label);
    fprintf('\n===== Measurement design %d/%d: %s =====\n',k,numel(designs),label);
    o=userOpts;
    o.measurement_label=label;
    o.measurement_wavenumbers_cm1=designs(k).wn;
    o.measurement_breath_times_h=designs(k).times;
    o.measurement_phase_s=designs(k).phases;
    R=run_reduced_patient_pinn([],fullfile(outputDir,char(label)),o);
    allResults{k}=R;

    row=add_design_columns(R.summary,R.parameters,label, ...
        numel(R.selected_wavenumbers_cm1),selected_count(designs(k).times,11), ...
        numel(R.selected_phase_indices));
    comparison=append_table(comparison,row);
    P=R.by_patient; P.design=repmat(label,height(P),1);
    P=movevars(P,'design','Before',1); byPatient=append_table(byPatient,P);
end

comparison.measurement_count=comparison.N_wavenumbers.* ...
    comparison.N_breath_times.*comparison.N_phase_samples;
comparison.measurement_reduction_percent=100*(1-comparison.measurement_count/ ...
    comparison.measurement_count(1));

validation=struct('summary',comparison,'by_patient',byPatient, ...
    'designs',designs,'results',{allResults},'user_options',userOpts);
writetable(comparison,fullfile(outputDir,'fisher_pinn_comparison.csv'));
writetable(byPatient,fullfile(outputDir,'fisher_pinn_by_patient.csv'));
save(fullfile(outputDir,'fisher_pinn_validation.mat'),'validation','-v7.3');
make_validation_plot(comparison,outputDir);
disp(comparison)
end

function n=selected_count(x,fallback)
if isempty(x); n=fallback; else; n=numel(x); end
end

function T=add_design_columns(S,P,label,nw,nt,np)
g=P.parameter=="Gamma_A"; gammaMAE=mean(abs(P.error(g)));
T=table(label,nw,nt,np,S.N_patients,S.N_temporal_holdout, ...
    S.BHB_RMSE_holdout_mM,S.BHB_MAE_holdout_mM, ...
    S.slope_RMSE_holdout_mM_h,S.direction_accuracy_holdout, ...
    S.sensitivity_rising,S.specificity_falling,S.balanced_accuracy_holdout, ...
    S.acetone_RMSE_holdout_ppm,gammaMAE, ...
    'VariableNames',{'design','N_wavenumbers','N_breath_times','N_phase_samples', ...
    'N_patients','N_temporal_holdout','BHB_RMSE_holdout_mM', ...
    'BHB_MAE_holdout_mM','slope_RMSE_holdout_mM_h', ...
    'direction_accuracy_holdout','sensitivity_rising','specificity_falling', ...
    'balanced_accuracy_holdout','acetone_RMSE_holdout_ppm','Gamma_A_MAE'});
end

function make_validation_plot(T,out)
f=figure('Color','w','Position',[100 100 1050 430]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile; bar(categorical(T.design),T.BHB_RMSE_holdout_mM); grid on
ylabel('BHB holdout RMSE (mmol/L)'); title('BHB reconstruction');
nexttile; bar(categorical(T.design),T.slope_RMSE_holdout_mM_h); grid on
ylabel('Slope holdout RMSE (mmol L^{-1} h^{-1})'); title('Trajectory reconstruction');
exportgraphics(f,fullfile(out,'fisher_pinn_comparison.png'),'Resolution',200);
close(f)
end

function T=append_table(T,A)
if isempty(T); T=A; else; T=[T;A]; end
end
