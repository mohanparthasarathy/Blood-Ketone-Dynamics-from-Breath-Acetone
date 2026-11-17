function U=run_fisher_uncertainty(fisherResults,dataFile,outputDir,userOpts)
%RUN_FISHER_UNCERTAINTY Profile Gamma_A under the selected Fisher D8 design.
% Research/synthetic-data workflow; alert labels are not a medical diagnosis.

here=fileparts(mfilename('fullpath')); projectDir=fileparts(here);
if nargin<1 || isempty(fisherResults)
    R=load(fullfile(projectDir,'fisher_pinn_validation','outputs', ...
        'fisher_validation_full','fisher_pinn_validation.mat'));
    labels=string(R.validation.summary.design);
    j=find(labels=="fisher_D8",1);
    if isempty(j); error('fisher_D8 results were not found.'); end
    fisherResults=R.validation.results{j};
end
if nargin<2 || isempty(dataFile)
    dataFile=fullfile(projectDir,'synthetic_data','outputs','cohort', ...
        'synthetic_cohort.mat');
end
if nargin<3 || isempty(outputDir); outputDir=fullfile(here,'outputs'); end
if nargin<4; userOpts=struct; end
o=options(userOpts); if ~exist(outputDir,'dir'); mkdir(outputDir); end
addpath(fullfile(projectDir,'fisher_pinn_validation'),'-begin');
addpath(fullfile(projectDir,'synthetic_data'),'-begin');
S=load(dataFile); ids=unique(fisherResults.predictions.patient_id);
if ~isempty(o.patient_ids); ids=ids(ismember(ids,o.patient_ids)); end
if isempty(ids); error('No requested patients found in PINN results.'); end

profile=[]; bands=[]; coverage=[];
profileFits=cell(numel(ids),1);
for ii=1:numel(ids)
    id=ids(ii); fprintf('\nGamma profile patient %d/%d (ID %d)\n',ii,numel(ids),id);
    center=parameter_value(fisherResults.parameters,id,"Gamma_A","estimate");
    truthGamma=parameter_value(fisherResults.parameters,id,"Gamma_A","truth");
    if isempty(o.gamma_grid_values)
        factors=linspace(1-o.gamma_span,1+o.gamma_span,o.gamma_grid_size);
        grid=unique(max(o.gamma_bounds(1),min(o.gamma_bounds(2),center*factors)));
    else
        grid=unique(max(o.gamma_bounds(1),min(o.gamma_bounds(2), ...
            double(o.gamma_grid_values(:).'))));
    end
    fits=cell(numel(grid),1); stat=zeros(numel(grid),1);
    for j=1:numel(grid)
        fprintf('  Gamma %.4f (%d/%d)\n',grid(j),j,numel(grid));
        po=design_options(fisherResults.options);
        po.gamma_mode="population"; po.gamma_population=grid(j); ...
        po.patient_ids=id; po.max_patients=1; ...
        po.max_iterations=o.profile_iterations; ...
        po.min_iterations=min(o.profile_min_iterations,o.profile_iterations); ...
        po.num_restarts=o.profile_restarts; ...
        po.refit_iterations=o.profile_refit_iterations; ...
        po.use_lbfgs=o.profile_use_lbfgs; ...
        po.lbfgs_iterations=o.profile_lbfgs_iterations; ...
        po.num_collocation=o.num_collocation; po.fast_nx=o.fast_nx; ...
        po.print_every=o.print_every; po.make_plots=false; ...
        po.save_results=false; po.random_seed=o.random_seed+1000*ii+j;
        cacheFile=fullfile(outputDir,sprintf( ...
            'profile_patient_%d_gamma_%.6f.mat',id,grid(j)));
        [useCache,C]=find_checkpoint(outputDir,id,grid(j),o.resume);
        if useCache
            fprintf('    loaded checkpoint\n'); fits{j}=C.fitResult;
        else
            fits{j}=run_reduced_patient_pinn(dataFile,fullfile(outputDir,'temporary'),po);
            fitResult=fits{j}; fitResult.trained_models={}; gamma_value=grid(j); %#ok<NASGU>
            save(cacheFile,'fitResult','gamma_value','-v7.3');
        end
        returnedIDs=unique(fits{j}.predictions.patient_id);
        if numel(returnedIDs)~=1 || returnedIDs~=id
            error(['Profile requested patient %d but run_patient_pinn returned ' ...
                'patient(s) %s. Check the standalone Fisher validation path.'], ...
                id,mat2str(returnedIDs(:).'));
        end
        P=fits{j}.predictions; m=P.used_for_fit;
        stat(j)=sum(((P.plateau_pred_ppm(m)-P.plateau_observed_ppm(m))/ ...
            o.acetone_sigma_ppm).^2);
    end
    delta=stat-min(stat); cutoff=2*gammaincinv(o.confidence_level,.5);
    accepted=delta<=cutoff; [~,bestj]=min(stat); profileFits{ii}=fits;
    lo=min(grid(accepted)); hi=max(grid(accepted));
    if accepted(1) || accepted(end)
        warning('Patient %d profile interval touches its Gamma grid boundary.',id);
    end
    scenario=string(fits{bestj}.predictions.scenario(1));
    T=table(repmat(id,numel(grid),1),repmat(scenario,numel(grid),1),grid(:), ...
        stat,delta,accepted,repmat(cutoff,numel(grid),1), ...
        'VariableNames',{'patient_id','scenario','Gamma_A','chi_square', ...
        'delta_chi_square','accepted','cutoff'});
    profile=append(profile,T);

    % Pointwise profile envelope at observed breath times.
    template=fits{bestj}.predictions; n=height(template);
    Bmat=nan(n,numel(grid)); Smat=Bmat;
    for j=1:numel(grid)
        Bmat(:,j)=fits{j}.predictions.BHB_pred_mM;
        Smat(:,j)=fits{j}.predictions.slope_pred_mM_h;
    end
    Blo=min(Bmat(:,accepted),[],2); Bhi=max(Bmat(:,accepted),[],2);
    Slo=min(Smat(:,accepted),[],2); Shi=max(Smat(:,accepted),[],2);
    status=strings(n,1);
    for k=1:n
        if Bhi(k)>=o.urgent_bhb_mM || ...
                (Bhi(k)>=o.monitor_bhb_mM && Shi(k)>=o.rapid_rise_mM_h)
            status(k)="urgent-evaluation";
        elseif Bhi(k)>=o.monitor_bhb_mM || Shi(k)>=o.rising_mM_h
            status(k)="monitor";
        else; status(k)="stable"; end
    end
    Bt=template.BHB_true_mM; St=template.slope_true_mM_h;
    TB=table(repmat(id,n,1),template.time_h,template.scenario,template.temporal_holdout, ...
        fits{bestj}.predictions.BHB_pred_mM,Blo,Bhi,Bt, ...
        fits{bestj}.predictions.slope_pred_mM_h,Slo,Shi,St,status, ...
        'VariableNames',{'patient_id','time_h','scenario','temporal_holdout', ...
        'BHB_point_mM','BHB_lower_mM','BHB_upper_mM','BHB_truth_mM', ...
        'slope_point_mM_h','slope_lower_mM_h','slope_upper_mM_h', ...
        'slope_truth_mM_h','stress_status'});
    bands=append(bands,TB);
    h=template.temporal_holdout;
    C=table(id,scenario,truthGamma,lo,hi,truthGamma>=lo && truthGamma<=hi, ...
        mean(Bt(h)>=Blo(h) & Bt(h)<=Bhi(h)), ...
        mean(St(h)>=Slo(h) & St(h)<=Shi(h)),mean(Bhi(h)-Blo(h)), ...
        'VariableNames',{'patient_id','scenario','Gamma_truth','Gamma_lower', ...
        'Gamma_upper','Gamma_covered','BHB_holdout_coverage', ...
        'slope_holdout_coverage','mean_BHB_band_width_mM'});
    coverage=append(coverage,C);
end

sensitivity=dimensionless_sensitivity(fisherResults,S,ids,o);
alerts=summarize_alerts(bands);
summary=table(o.confidence_level,height(coverage),mean(coverage.Gamma_covered), ...
    mean(coverage.BHB_holdout_coverage),mean(coverage.slope_holdout_coverage), ...
    mean(coverage.mean_BHB_band_width_mM), ...
    'VariableNames',{'confidence_level','N_patients','Gamma_coverage', ...
    'BHB_holdout_coverage','slope_holdout_coverage','mean_BHB_band_width_mM'});
U=struct('design',"fisher_D8",'summary',summary,'coverage',coverage,'profile',profile,'bands',bands, ...
    'alerts',alerts,'sensitivity',sensitivity,'options',o,'profile_fits',{profileFits});
writetable(summary,fullfile(outputDir,'uncertainty_summary.csv'));
writetable(coverage,fullfile(outputDir,'uncertainty_coverage.csv'));
writetable(profile,fullfile(outputDir,'gamma_profiles.csv'));
writetable(bands,fullfile(outputDir,'state_uncertainty_bands.csv'));
writetable(sensitivity,fullfile(outputDir,'dimensionless_sensitivity.csv'));
writetable(alerts,fullfile(outputDir,'upper_confidence_alerts.csv'));
save(fullfile(outputDir,'pinn_uncertainty_results.mat'),'U','-v7.3');
make_plots(U,outputDir); disp(summary); disp(coverage)
end

function [found,C]=find_checkpoint(outputDir,id,gamma,resume)
found=false; C=struct;
if ~resume; return; end
files=dir(fullfile(outputDir,sprintf('profile_patient_%d_*.mat',id)));
for k=1:numel(files)
    candidate=load(fullfile(files(k).folder,files(k).name));
    if isfield(candidate,'gamma_value') && isfield(candidate,'fitResult') && ...
            abs(candidate.gamma_value-gamma)<1e-9
        found=true; C=candidate; return
    end
end
end

function A=summarize_alerts(B)
A=[]; ids=unique(B.patient_id);
for id=ids.'
    Q=B(B.patient_id==id,:); monitor=Q.stress_status=="monitor" | ...
        Q.stress_status=="urgent-evaluation";
    urgent=Q.stress_status=="urgent-evaluation";
    firstMonitor=NaN; firstUrgent=NaN;
    if any(monitor); firstMonitor=min(Q.time_h(monitor)); end
    if any(urgent); firstUrgent=min(Q.time_h(urgent)); end
    row=table(id,string(Q.scenario(1)),firstMonitor,firstUrgent, ...
        max(Q.BHB_upper_mM),max(Q.slope_upper_mM_h),sum(monitor),sum(urgent), ...
        'VariableNames',{'patient_id','scenario','first_monitor_time_h', ...
        'first_urgent_time_h','maximum_BHB_upper_mM', ...
        'maximum_slope_upper_mM_h','N_monitor_or_higher','N_urgent'});
    A=append(A,row);
end
end

function T=dimensionless_sensitivity(R,S,ids,o)
T=[]; factors=o.sensitivity_factors(:); base=S.base;
for ii=1:numel(ids)
    id=ids(ii); M=S.meta(S.meta.patient_id==id,:); [~,ix]=sort(M.breath_time_h); M=M(ix,:);
    P=R.predictions(R.predictions.patient_id==id,:); [~,ix]=sort(P.time_h); P=P(ix,:);
    gamma=parameter_value(R.parameters,id,"Gamma_A","estimate");
    phaseTimes=selected_phase_times(S,M,R.options.measurement_phase_s);
    a=P.plateau_pred_ppm./max(eps,pulmonary_gain_vector(M,base,o.fast_nx,phaseTimes))./ ...
        (base.slow.kappa_QA*base.scales.K_star_mM/gamma/base.scales.lambda_ba*24450);
    baseline=predicted_plateau(M,base,a,gamma,1,1,1,o.fast_nx,phaseTimes);
    for name=["Pe","Da_aw","Gamma_A"]
        for j=1:numel(factors)
            fp=1; fd=1; fg=1;
            if name=="Pe"; fp=factors(j); elseif name=="Da_aw"; fd=factors(j); else; fg=factors(j); end
            y=predicted_plateau(M,base,a,gamma,fp,fd,fg,o.fast_nx,phaseTimes);
            row=table(id,string(M.scenario(1)),name,factors(j),mean(y), ...
                100*mean((y-baseline)./max(abs(baseline),eps)), ...
                'VariableNames',{'patient_id','scenario','group','factor', ...
                'mean_plateau_ppm','mean_percent_change'});
            T=append(T,row);
        end
    end
end
end

function y=predicted_plateau(M,base,a,gamma,fp,fd,fg,nx,phaseTimes)
gain=zeros(height(M),1); p=base; p.fast.nx=nx;
for k=1:height(M)
    p.fast.Pe=M.Pe(k)*fp; p.fast.Da_aw=M.Da_aw(k)*fd;
    z=solve_fast_breath(1,p,M.U_scale(k));
    gain(k)=mean(interp1(z.time_s,z.c(end,:),phaseTimes(:,k),'pchip'));
end
g=gamma*fg; cstar=base.slow.kappa_QA*base.scales.K_star_mM/g/ ...
    base.scales.lambda_ba*24450; y=gain.*cstar.*a;
end
function g=pulmonary_gain_vector(M,base,nx,phaseTimes)
g=zeros(height(M),1); p=base; p.fast.nx=nx;
for k=1:height(M); p.fast.Pe=M.Pe(k); p.fast.Da_aw=M.Da_aw(k);
 z=solve_fast_breath(1,p,M.U_scale(k));
 g(k)=mean(interp1(z.time_s,z.c(end,:),phaseTimes(:,k),'pchip')); end
end

function T=selected_phase_times(S,M,requested)
allTimes=double(S.time_s(:,1));
if isempty(requested); idx=(1:numel(allTimes)).';
else
    idx=zeros(numel(requested),1);
    for k=1:numel(requested); [~,idx(k)]=min(abs(allTimes-requested(k))); end
    idx=unique(idx,'stable');
end
rows=zeros(height(M),1);
for k=1:height(M)
    candidates=find(S.meta.patient_id==M.patient_id(k));
    [~,j]=min(abs(S.meta.breath_time_h(candidates)-M.breath_time_h(k)));
    rows(k)=candidates(j);
end
T=double(S.time_s(idx,rows));
end

function po=design_options(base)
fields={'measurement_label','measurement_wavenumbers_cm1', ...
    'measurement_phase_s','measurement_breath_times_h', ...
    'fit_spectral_nuisance','add_spectral_nuisance', ...
    'nuisance_scale_abs','nuisance_random_seed','gamma_prior_cv'};
po=struct;
for k=1:numel(fields)
    if isfield(base,fields{k}); po.(fields{k})=base.(fields{k}); end
end
end
function x=parameter_value(T,id,name,column)
row=T.patient_id==id & string(T.parameter)==name; x=T.(column)(row);
if isempty(x); error('Parameter %s not found for patient %d.',name,id); end
x=x(1);
end
function T=append(T,A); if isempty(T); T=A; else; T=[T;A]; end; end

function make_plots(U,out)
ids=unique(U.bands.patient_id); f=figure('Color','w','Position',[80 80 1200 650]);
tiledlayout(2,ceil(numel(ids)/2),'TileSpacing','compact');
for k=1:numel(ids); nexttile; Q=U.bands(U.bands.patient_id==ids(k),:);
 fill([Q.time_h;flipud(Q.time_h)],[Q.BHB_lower_mM;flipud(Q.BHB_upper_mM)], ...
     [.75 .87 1],'EdgeColor','none'); hold on
 plot(Q.time_h,Q.BHB_truth_mM,'k-','LineWidth',1.5);
 plot(Q.time_h,Q.BHB_point_mM,'b-','LineWidth',1.5); grid on
 title(sprintf('Patient %d: %s',ids(k),string(Q.scenario(1)))); xlabel('Time (h)'); ylabel('BHB (mmol/L)');
end
exportgraphics(f,fullfile(out,'BHB_uncertainty_bands.png'),'Resolution',200); close(f)
end

function o=options(o)
d=struct('patient_ids',[],'confidence_level',.95,'gamma_span',.45, ...
 'gamma_grid_values',[], ...
 'gamma_grid_size',13,'gamma_bounds',[.3 4.5],'acetone_sigma_ppm',.04, ...
 'profile_iterations',5000,'profile_min_iterations',2500,'profile_restarts',1, ...
 'profile_refit_iterations',500,'profile_use_lbfgs',true, ...
 'profile_lbfgs_iterations',100,'num_collocation',121,'fast_nx',81, ...
 'print_every',1000,'random_seed',73001,'sensitivity_factors',[.5 .75 1 1.25 1.5], ...
 'monitor_bhb_mM',1.5,'urgent_bhb_mM',3.0,'rising_mM_h',.15, ...
 'rapid_rise_mM_h',.30,'resume',true);
f=fieldnames(d); for k=1:numel(f); if ~isfield(o,f{k}); o.(f{k})=d.(f{k}); end; end
end
