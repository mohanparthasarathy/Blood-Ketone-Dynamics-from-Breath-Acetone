function results = run_patient_pinn(dataFile,outputDir,userOpts)
%RUN_PATIENT_PINN Train one metabolic PINN per held-out virtual patient.
% Requires MATLAB R2024b+ and Deep Learning Toolbox (dljacobian).

here=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(dataFile)
    projectDir=fileparts(here);
    dataFile=fullfile(projectDir,'synthetic_data','outputs','cohort', ...
        'synthetic_cohort.mat');
end
if nargin<2 || isempty(outputDir); outputDir=fullfile(here,'outputs'); end
if nargin<3; userOpts=struct; end
o=pinn_options(userOpts);
if ~isfile(dataFile); error('Cohort file not found: %s',dataFile); end
if ~exist(outputDir,'dir'); mkdir(outputDir); end
if exist('dljacobian','file')~=2
    error('dljacobian requires MATLAB R2024b or newer and Deep Learning Toolbox.');
end

syntheticDir=fullfile(fileparts(here),'synthetic_data');
if isfolder(syntheticDir); addpath(syntheticDir); end
if exist('solve_fast_breath','file')~=2
    error('Add the sibling synthetic_data folder to the MATLAB path.');
end
S=load(dataFile);
if ~isfield(S,'cohort_version') || S.cohort_version<2
    error('Regenerate the corrected Version 2 synthetic cohort first.');
end

decoded=decode_pinn_spectra(S.spectra,S.cleanSpectra,S.wn,S.base);
nTail=max(2,round(.1*size(decoded.C_hat_ppm,1)));
CobsAll=mean(decoded.C_hat_ppm(end-nTail:end,:),1).';
CtrueAll=mean(decoded.C_true_ppm(end-nTail:end,:),1).';

ids = unique(S.meta.patient_id(string(S.meta.split)=="test"));

if ~isempty(o.patient_ids)
    ids = ids(ismember(ids,o.patient_ids));

    if isempty(ids)
        error('None of the requested patient_ids are in the test split.');
    end
end

if isfinite(o.max_patients)
    ids = ids(1:min(numel(ids),o.max_patients));
end
predictions=[]; parameters=[]; patients=[]; trained=cell(numel(ids),1);
for ip=1:numel(ids)
    id=ids(ip); idx=find(S.meta.patient_id==id);
    [~,ord]=sort(S.meta.breath_time_h(idx)); idx=idx(ord); M=S.meta(idx,:);
    fprintf('\nPINN patient %d/%d (ID %d, %s)\n',ip,numel(ids),id,string(M.scenario(1)));
    gain=dimensionless_gains(M,S.base,o.fast_nx);
    t=M.breath_time_h; Cobs=CobsAll(idx);
    fitMask=mod((1:numel(t)).',o.holdout_every)~=0;
    truth=true_patient_parameters(M);
    fit=train_patient_pinn(t,Cobs,gain,S.base,fitMask,o,truth,string(M.scenario(1)));
    trained{ip}=fit;

    [state,dstate]=evaluate_pinn(fit.net,t,o.t_max_h);
    B=state(1,:).'; slope=dstate(1,:).';
    Cstar=S.base.slow.kappa_QA*S.base.scales.K_star_mM/fit.Gamma_A ...
        /S.base.scales.lambda_ba*24450;
    Cpred=gain.*Cstar.*state(3,:).'; holdout=~fitMask;
    P=table(repmat(id,numel(t),1),t,string(M.scenario),fitMask,holdout, ...
        M.BHB_mM,B,M.dBHB_dt_mM_h,slope,Cobs,Cpred,CtrueAll(idx), ...
        'VariableNames',{'patient_id','time_h','scenario','used_for_fit', ...
        'temporal_holdout','BHB_true_mM','BHB_pred_mM','slope_true_mM_h', ...
        'slope_pred_mM_h','plateau_observed_ppm','plateau_pred_ppm', ...
        'plateau_clean_ppm'});
    predictions=append_table(predictions,P);

    names=["r_pre","r_post","tau_c","k_r","Gamma_A"];
    est=[fit.r_pre fit.r_post fit.tau_c fit.k_r fit.Gamma_A];
    target=truth(4:8);
    T=table(repmat(id,5,1),names.',est(:),target(:),est(:)-target(:), ...
        'VariableNames',{'patient_id','parameter','estimate','truth','error'});
    parameters=append_table(parameters,T);
    h=holdout;
    row=table(id,string(M.scenario(1)),sum(fitMask),sum(h),fit.final_loss, ...
        sqrt(mean((B(h)-M.BHB_mM(h)).^2)), ...
        sqrt(mean((slope(h)-M.dBHB_dt_mM_h(h)).^2)), ...
        sqrt(mean((Cpred(h)-Cobs(h)).^2)), ...
        'VariableNames',{'patient_id','scenario','N_fit','N_holdout','final_loss', ...
        'BHB_RMSE_holdout_mM','slope_RMSE_holdout_mM_h','acetone_RMSE_holdout_ppm'});
    patients=append_table(patients,row);
end

summary=summarize_pinn(predictions,patients,o.gamma_mode);
results = struct( ...
    'summary',summary, ...
    'by_patient',patients, ...
    'parameters',parameters, ...
    'predictions',predictions, ...
    'trained_models',{trained}, ...
    'options',o);

if o.save_results
    writetable(predictions, ...
        fullfile(outputDir,'pinn_predictions.csv'));

    writetable(parameters, ...
        fullfile(outputDir,'pinn_parameters.csv'));

    writetable(patients, ...
        fullfile(outputDir,'pinn_by_patient.csv'));

    writetable(summary, ...
        fullfile(outputDir,'pinn_summary.csv'));

    save(fullfile(outputDir,'patient_pinn_results.mat'), ...
        'results','-v7.3');
end

if o.make_plots
    make_pinn_plots(predictions,outputDir,o.gamma_mode);
end

disp(summary)
end

function gain=dimensionless_gains(M,base,nx)
p=base; p.fast.nx=nx; gain=zeros(height(M),1);
for j=1:height(M)
    p.fast.Pe=M.Pe(j); p.fast.Da_aw=M.Da_aw(j);
    unit=solve_fast_breath(1,p,M.U_scale(j));
    n=max(2,round(.1*numel(unit.time_s)));
    gain(j)=mean(unit.c(end,end-n:end));
end
end

function truth=true_patient_parameters(M)
s=string(M.scenario(1)); r0=M.r0(1); dr=M.delta_r(1);
if s=="developing"; rp=r0; rq=r0+dr;
elseif s=="resolving"; rp=r0+dr; rq=r0;
else; rp=r0; rq=r0; end
truth=[.18 .06 .06 rp rq M.tau_c(1) M.k_r(1) M.Gamma_A(1)];
end

function S=summarize_pinn(P,byPatient,mode)
h=P.temporal_holdout; changing=h & string(P.scenario)~="baseline";
truth=P.slope_true_mM_h(changing)>0; guess=P.slope_pred_mM_h(changing)>0;
sens=safe_mean(guess(truth)); spec=safe_mean(~guess(~truth));
S=table(string(mode),numel(unique(P.patient_id)),sum(h), ...
    sqrt(mean((P.BHB_pred_mM(h)-P.BHB_true_mM(h)).^2)), ...
    mean(abs(P.BHB_pred_mM(h)-P.BHB_true_mM(h))), ...
    sqrt(mean((P.slope_pred_mM_h(h)-P.slope_true_mM_h(h)).^2)), ...
    mean(guess==truth),sens,spec,mean([sens spec],'omitnan'), ...
    mean(byPatient.acetone_RMSE_holdout_ppm), ...
    'VariableNames',{'gamma_mode','N_patients','N_temporal_holdout', ...
    'BHB_RMSE_holdout_mM','BHB_MAE_holdout_mM','slope_RMSE_holdout_mM_h', ...
    'direction_accuracy_holdout','sensitivity_rising','specificity_falling', ...
    'balanced_accuracy_holdout','acetone_RMSE_holdout_ppm'});
end
function T=append_table(T,A); if isempty(T); T=A; else; T=[T;A]; end; end
function m=safe_mean(x); if isempty(x); m=NaN; else; m=mean(x); end; end

function o=pinn_options(o)
userSetIterations=isfield(o,'iterations');
userSetMaxIterations=isfield(o,'max_iterations');
d=struct('patient_ids',[], ...
'save_results',true, ...
'make_plots',true, ...
'gamma_mode',"oracle",'gamma_population',1.5,'gamma_prior_cv',.10, ...
 'max_patients',Inf,'holdout_every',3,'fast_nx',121,'t_max_h',8, ...
 'num_collocation',161,'hidden_width',48,'hidden_layers',3, ...
 'iterations',4000,'max_iterations',8000,'min_iterations',3000, ...
 'learning_rate',1e-3,'print_every',500,'num_restarts',3, ...
 'refit_iterations',1000,'curriculum_data_fraction',.15, ...
 'curriculum_ramp_fraction',.35,'early_stop_check_every',250, ...
 'early_stop_patience',1250,'early_stop_rel_tol',2e-3, ...
 'lambda_data',10,'lambda_physics',1,'lambda_ic',2, ...
 'lambda_positive',1,'lambda_gamma',.002,'lambda_forcing',2e-3, ...
 'r_low_prior_center',.20,'r_high_prior_center',2.5,'r_prior_scale',1.5, ...
 'kr_prior_center',2,'kr_prior_scale',1.5, ...
 'ic_target',[.18 .06 .06],'selection_physics_weight',.1, ...
 'use_lbfgs',true,'lbfgs_iterations',300,'lbfgs_print_every',50, ...
 'lbfgs_gradient_tolerance',1e-5,'lbfgs_step_tolerance',1e-6, ...
 'verbose',true,'random_seed',24680);
f=fieldnames(d); for k=1:numel(f); if ~isfield(o,f{k}); o.(f{k})=d.(f{k}); end; end
if userSetIterations && ~userSetMaxIterations; o.max_iterations=o.iterations; end
end
