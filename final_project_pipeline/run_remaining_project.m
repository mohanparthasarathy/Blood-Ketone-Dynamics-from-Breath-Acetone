function final = run_remaining_project(userOpts)
%RUN_REMAINING_PROJECT Resumable final experiments for the breath-ketone project.
%
% final = run_remaining_project
% final = run_remaining_project(struct('mode',"smoke"))
%
% The runner never modifies the established patient_pinn or Fisher folders.
% Every expensive condition is saved independently and automatically reused.
% Alert labels are research outputs and are not medical diagnoses.

if nargin<1; userOpts=struct; end
here=fileparts(mfilename('fullpath')); root=fileparts(here);
o=defaults(userOpts,here);
ensure_folder(o.output_dir);

required={fullfile(root,'synthetic_data'), ...
    fullfile(root,'fisher_pinn_validation')};
for k=1:numel(required)
    if ~isfolder(required{k}); error('Required sibling folder missing: %s',required{k}); end
    addpath(required{k},'-begin');
end

dataFile=resolve_data_file(root,o.data_file);
fprintf('\nFINAL PROJECT PIPELINE\nOutput: %s\n',o.output_dir);
fprintf('Completed condition folders will be reused automatically.\n');

% The experiment list changes one feature at a time around the reference.
conditions=make_conditions(o);
designs=make_designs();
manifest=[]; allPatient=[]; allPredictions=[]; allParameters=[];

for c=1:numel(conditions)
    C=conditions(c);
    conditionData=prepare_noise_data(dataFile,o.output_dir,C,o);
    for d=1:numel(designs)
        D=designs(d);
        tag=sprintf('%s__%s',char(C.label),char(D.label));
        runDir=fullfile(o.output_dir,'checkpoints',tag);
        resultFile=fullfile(runDir,'patient_pinn_results.mat');
        ensure_folder(runDir);

        if isfile(resultFile) && o.resume
            fprintf('\n[%s] loading checkpoint\n',tag);
            Z=load(resultFile); R=Z.results;
        else
            fprintf('\n[%s] fitting\n',tag);
            po=pinn_options(o,C,D);
            R=run_reduced_patient_pinn(conditionData,runDir,po);
        end

        row=summary_row(R,C,D);
        manifest=append(manifest,row);
        P=R.by_patient; P.condition=repmat(C.label,height(P),1);
        P.design=repmat(D.label,height(P),1);
        P=movevars(P,{'condition','design'},'Before',1);
        allPatient=append(allPatient,P);
        P=R.predictions; P.condition=repmat(C.label,height(P),1);
        P.design=repmat(D.label,height(P),1);
        P=movevars(P,{'condition','design'},'Before',1);
        allPredictions=append(allPredictions,P);
        P=R.parameters; P.condition=repmat(C.label,height(P),1);
        P.design=repmat(D.label,height(P),1);
        P=movevars(P,{'condition','design'},'Before',1);
        allParameters=append(allParameters,P);

        partial=struct('summary',manifest,'by_patient',allPatient, ...
            'predictions',allPredictions,'parameters',allParameters, ...
            'options',o,'completed_at',datetime('now'));
        save(fullfile(o.output_dir,'pipeline_checkpoint.mat'),'partial','-v7.3');
        writetable(manifest,fullfile(o.output_dir,'progress_summary.csv'));
    end
end

alerts=evaluate_alerts(allPredictions,o);
regimes=dimensionless_regime_map(dataFile,o);
uncertaintyAlerts=load_uncertainty_alerts(root,o);
conclusions=make_conclusions(manifest,alerts.summary,regimes);

final=struct('summary',manifest,'by_patient',allPatient, ...
    'predictions',allPredictions,'parameters',allParameters, ...
    'alert_summary',alerts.summary,'alert_events',alerts.events, ...
    'regime_map',regimes,'uncertainty_alerts',uncertaintyAlerts, ...
    'conclusions',conclusions,'conditions',conditions,'designs',designs, ...
    'options',o,'completed_at',datetime('now'));

writetable(manifest,fullfile(o.output_dir,'robustness_summary.csv'));
writetable(allPatient,fullfile(o.output_dir,'robustness_by_patient.csv'));
writetable(alerts.summary,fullfile(o.output_dir,'alert_validation_summary.csv'));
writetable(alerts.events,fullfile(o.output_dir,'alert_patient_events.csv'));
writetable(regimes,fullfile(o.output_dir,'dimensionless_regime_map.csv'));
writetable(conclusions,fullfile(o.output_dir,'final_conclusions.csv'));
save(fullfile(o.output_dir,'final_project_results.mat'),'final','-v7.3');
make_figures(final,o.output_dir);

fprintf('\nFINAL PIPELINE COMPLETE\n');
disp(final.summary)
disp(final.alert_summary)
fprintf('Saved to %s\n',o.output_dir);
end

function o=defaults(o,here)
d=struct('mode',"full",'output_dir',fullfile(here,'outputs','final_run'), ...
 'data_file',"",'resume',true,'random_seed',88001, ...
 'noise_multipliers',[.5 1 2],'breath_counts',[4 6 8], ...
 'phase_counts',[8 15],'max_patients',Inf, ...
 'max_iterations',8000,'min_iterations',3000,'num_restarts',3, ...
 'refit_iterations',1000,'use_lbfgs',true,'lbfgs_iterations',300, ...
 'print_every',1000,'gamma_mode',"prior",'gamma_prior_cv',.10, ...
 'monitor_bhb_mM',1.5,'urgent_bhb_mM',3.0, ...
 'rising_mM_h',.15,'rapid_rise_mM_h',.30, ...
 'Pe_values',[12 20 30 45 60],'Da_values',[.35 .6 .9 1.2 1.5], ...
 'Gamma_values',[.7 1.2 1.8 2.4 2.8], ...
 'uncertainty_file',"");
f=fieldnames(d); for k=1:numel(f); if ~isfield(o,f{k}); o.(f{k})=d.(f{k}); end; end
if string(o.mode)=="smoke"
    o.noise_multipliers=[1 2]; o.breath_counts=[4 8]; o.phase_counts=8;
    o.max_patients=1; o.max_iterations=1500; o.min_iterations=500;
    o.num_restarts=1; o.refit_iterations=100; o.use_lbfgs=false;
    o.Pe_values=[12 30 60]; o.Da_values=[.35 .9 1.5];
    o.Gamma_values=[.7 1.8 2.8];
end
end

function file=resolve_data_file(root,requested)
if strlength(string(requested))>0; file=char(requested); else
    file=fullfile(root,'synthetic_data','outputs','cohort','synthetic_cohort.mat');
end
if ~isfile(file); error('Synthetic cohort not found: %s',file); end
end

function C=make_conditions(o)
C=struct('label',{},'noise_multiplier',{},'breath_count',{},'phase_count',{});
C(end+1)=entry("reference",1,8,15);
for x=o.noise_multipliers
    if x~=1; C(end+1)=entry("noise_"+safe_num(x),x,8,15); end
end
for x=o.breath_counts
    if x~=8; C(end+1)=entry("breaths_"+string(x),1,x,15); end
end
for x=o.phase_counts
    if x~=15; C(end+1)=entry("phases_"+string(x),1,8,x); end
end
end
function x=entry(label,n,b,p)
x=struct('label',string(label),'noise_multiplier',n,'breath_count',b,'phase_count',p);
end
function s=safe_num(x); s=replace(string(x),".","p"); end

function D=make_designs()
D=struct( ...
 'label',{"fisher_D8","uniform_8"}, ...
 'wn',{[1210 1226 1450 1630 1720 1736 1750 1900], ...
       [1000 1136 1256 1390 1510 1646 1766 1900]});
end

function file=prepare_noise_data(source,out,C,o)
if C.noise_multiplier==1; file=source; return; end
noiseDir=fullfile(out,'datasets'); ensure_folder(noiseDir);
file=fullfile(noiseDir,char(C.label+".mat"));
if isfile(file) && o.resume; return; end
S=load(source);
rng(o.random_seed+round(100*C.noise_multiplier),'twister');
sigma=reshape(double(S.meta.noise_sd_abs)*C.noise_multiplier,1,1,[]);
S.spectra=single(double(S.cleanSpectra)+randn(size(S.cleanSpectra)).*sigma);
S.meta.noise_sd_abs=S.meta.noise_sd_abs*C.noise_multiplier;
spectra=S.spectra; cleanSpectra=S.cleanSpectra; time_s=S.time_s; wn=S.wn;
meta=S.meta; opts=S.opts; base=S.base; cohort_version=S.cohort_version; %#ok<NASGU>
save(file,'spectra','cleanSpectra','time_s','wn','meta','opts','base', ...
    'cohort_version','-v7.3');
end

function po=pinn_options(o,C,D)
baseTimes=[.5 1.25 2 2.75 3.5 5 5.75 8];
basePhases=[1.35 1.85 2.50 3.15 3.35 3.50 3.65 3.85 ...
    4.00 4.15 4.35 4.50 4.65 4.85 5.00];
ti=unique(round(linspace(1,numel(baseTimes),C.breath_count)));
pi=unique(round(linspace(1,numel(basePhases),C.phase_count)));
po=struct('measurement_label',C.label+"__"+string(D.label), ...
 'measurement_wavenumbers_cm1',D.wn, ...
 'measurement_breath_times_h',baseTimes(ti), ...
 'measurement_phase_s',basePhases(pi), ...
 'gamma_mode',o.gamma_mode,'gamma_prior_cv',o.gamma_prior_cv, ...
 'max_patients',o.max_patients,'max_iterations',o.max_iterations, ...
 'min_iterations',o.min_iterations,'num_restarts',o.num_restarts, ...
 'refit_iterations',o.refit_iterations,'use_lbfgs',o.use_lbfgs, ...
 'lbfgs_iterations',o.lbfgs_iterations,'print_every',o.print_every, ...
 'make_plots',false,'save_results',true, ...
 'random_seed',o.random_seed+sum(double(char(C.label)))+sum(double(char(D.label))));
end

function T=summary_row(R,C,D)
S=R.summary;
T=table(C.label,string(D.label),C.noise_multiplier,C.breath_count,C.phase_count, ...
 S.N_patients,S.N_temporal_holdout,S.BHB_RMSE_holdout_mM,S.BHB_MAE_holdout_mM, ...
 S.slope_RMSE_holdout_mM_h,S.direction_accuracy_holdout, ...
 S.balanced_accuracy_holdout,S.acetone_RMSE_holdout_ppm, ...
 'VariableNames',{'condition','design','noise_multiplier','N_breaths', ...
 'N_phases','N_patients','N_holdout','BHB_RMSE_mM','BHB_MAE_mM', ...
 'slope_RMSE_mM_h','direction_accuracy','balanced_accuracy','acetone_RMSE_ppm'});
end

function A=evaluate_alerts(P,o)
events=[]; groups=unique(P(:,{'condition','design'}),'rows');
for g=1:height(groups)
    m=P.condition==groups.condition(g) & P.design==groups.design(g);
    ids=unique(P.patient_id(m));
    for id=ids.'
        Q=sortrows(P(m & P.patient_id==id,:),'time_h');
        use=Q.temporal_holdout;
        trueUrg=urgent(Q.BHB_true_mM,Q.slope_true_mM_h,o);
        predUrg=urgent(Q.BHB_pred_mM,Q.slope_pred_mM_h,o);
        trueMon=monitor(Q.BHB_true_mM,Q.slope_true_mM_h,o)|trueUrg;
        predMon=monitor(Q.BHB_pred_mM,Q.slope_pred_mM_h,o)|predUrg;
        row=table(groups.condition(g),groups.design(g),id,string(Q.scenario(1)), ...
            first_time(Q.time_h,trueMon),first_time(Q.time_h,predMon), ...
            first_time(Q.time_h,trueUrg),first_time(Q.time_h,predUrg), ...
            mean(predUrg(use) & ~trueUrg(use)),mean(~predUrg(use) & trueUrg(use)), ...
            mean(predUrg(use)==trueUrg(use)),last_time(Q.time_h,predMon), ...
            'VariableNames',{'condition','design','patient_id','scenario', ...
            'truth_first_monitor_h','pred_first_monitor_h','truth_first_urgent_h', ...
            'pred_first_urgent_h','false_urgent_rate','missed_urgent_rate', ...
            'urgent_accuracy','pred_last_monitor_or_higher_h'});
        events=append(events,row);
    end
end
summary=[];
for g=1:height(groups)
    Q=events(events.condition==groups.condition(g)&events.design==groups.design(g),:);
    base=Q.scenario=="baseline"; dev=Q.scenario=="developing"; res=Q.scenario=="resolving";
    lead=Q.pred_first_urgent_h(dev)-Q.truth_first_urgent_h(dev);
    row=table(groups.condition(g),groups.design(g),mean(Q.false_urgent_rate), ...
        mean(Q.missed_urgent_rate),mean(Q.urgent_accuracy), ...
        mean(Q.false_urgent_rate(base),'omitnan'),mean(lead,'omitnan'), ...
        mean(Q.pred_last_monitor_or_higher_h(res),'omitnan'), ...
        'VariableNames',{'condition','design','false_urgent_rate', ...
        'missed_urgent_rate','urgent_accuracy','baseline_false_urgent_rate', ...
        'developing_detection_delay_h','resolving_last_alert_h'});
    summary=append(summary,row);
end
A=struct('summary',summary,'events',events);
end
function y=urgent(B,S,o); y=B>=o.urgent_bhb_mM | (B>=o.monitor_bhb_mM & S>=o.rapid_rise_mM_h); end
function y=monitor(B,S,o); y=B>=o.monitor_bhb_mM | S>=o.rising_mM_h; end
function x=first_time(t,m); if any(m); x=min(t(m)); else; x=NaN; end; end
function x=last_time(t,m); if any(m); x=max(t(m)); else; x=NaN; end; end

function T=dimensionless_regime_map(dataFile,o)
S=load(dataFile,'base'); base=S.base; T=[];
for G=o.Gamma_values
 for Pe=o.Pe_values
  for Da=o.Da_values
   p=base; p.slow.Gamma_A=G; p.fast.Pe=Pe; p.fast.Da_aw=Da;
   p.scales.A_star_mM=p.slow.kappa_QA*p.scales.K_star_mM/G;
   p.scales.C_star_ppm=p.scales.A_star_mM/p.scales.lambda_ba*24450;
   b=solve_fast_breath(1,p,1); c=b.C_exh_ppm;
   plateau=mean(c(end-max(2,round(.1*numel(c))):end));
   row=table(Pe,Da,G,plateau,trapz(b.time_s,c),max(c), ...
    'VariableNames',{'Pe','Da_aw','Gamma_A','plateau_ppm', ...
    'AUC_ppm_s','maximum_ppm'});
   T=append(T,row);
  end
 end
end
end

function T=load_uncertainty_alerts(root,o)
T=table;
if strlength(string(o.uncertainty_file))>0; file=char(o.uncertainty_file); else
 candidates={fullfile(root,'fisher_uncertainty','outputs','fisher_uncertainty_full', ...
     'pinn_uncertainty_results.mat'), ...
     fullfile(root,'fisher_uncertainty','outputs','pinn_uncertainty_results.mat')};
 file=''; for k=1:numel(candidates); if isfile(candidates{k}); file=candidates{k}; break; end; end
end
if isempty(file) || ~isfile(file)
    warning('Uncertainty result not found; uncertainty-alert table was skipped.'); return
end
Z=load(file); names=fieldnames(Z); U=[];
for k=1:numel(names); if isstruct(Z.(names{k})) && isfield(Z.(names{k}),'alerts'); U=Z.(names{k}); break; end; end
if ~isempty(U); T=U.alerts; end
end

function T=make_conclusions(S,A,R)
ref=S(S.condition=="reference",:); T=table;
for d=unique(S.design).'
 q=ref(ref.design==d,:); if isempty(q); continue; end
 row=table(d,q.BHB_RMSE_mM,q.slope_RMSE_mM_h,q.acetone_RMSE_ppm, ...
    A.baseline_false_urgent_rate(A.condition=="reference"&A.design==d), ...
    max(R.plateau_ppm)/min(R.plateau_ppm), ...
    'VariableNames',{'design','reference_BHB_RMSE_mM', ...
    'reference_slope_RMSE_mM_h','reference_acetone_RMSE_ppm', ...
    'baseline_false_urgent_rate','regime_plateau_max_to_min'});
 T=append(T,row);
end
end

function make_figures(F,out)
f=figure('Color','w','Position',[100 100 1100 430]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile; hold on
for d=unique(F.summary.design).'
 q=F.summary(F.summary.design==d,:);
 plot(categorical(q.condition),q.BHB_RMSE_mM,'-o','DisplayName',d);
end
grid on; ylabel('BHB holdout RMSE (mmol/L)'); title('Robustness'); legend('Location','best');
nexttile; scatter3(F.regime_map.Pe,F.regime_map.Da_aw,F.regime_map.plateau_ppm, ...
 45,F.regime_map.Gamma_A,'filled'); grid on; colorbar
xlabel('Pe'); ylabel('Da_{aw}'); zlabel('Plateau (ppm)'); title('Dimensionless regimes');
exportgraphics(f,fullfile(out,'final_robustness_and_regimes.png'),'Resolution',220);
close(f)
end

function ensure_folder(x); if ~exist(x,'dir'); mkdir(x); end; end
function T=append(T,A); if isempty(T); T=A; else; T=[T;A]; end; end
