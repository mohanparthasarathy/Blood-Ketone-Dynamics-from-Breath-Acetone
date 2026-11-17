function U=run_pinn_uncertainty(pinnResults,dataFile,outputDir,userOpts)
%RUN_PINN_UNCERTAINTY Profile Gamma_A, build state bands, sensitivities, alerts.
% Research/synthetic-data workflow; alert labels are not a medical diagnosis.

here=fileparts(mfilename('fullpath')); projectDir=fileparts(here);
if nargin<1 || isempty(pinnResults)
    R=load(fullfile(projectDir,'patient_pinn','outputs','prior10_v2', ...
        'patient_pinn_results.mat')); pinnResults=R.results;
end
if nargin<2 || isempty(dataFile)
    dataFile=fullfile(projectDir,'synthetic_data','outputs','cohort', ...
        'synthetic_cohort.mat');
end
if nargin<3 || isempty(outputDir); outputDir=fullfile(here,'outputs'); end
if nargin<4; userOpts=struct; end
o=options(userOpts); if ~exist(outputDir,'dir'); mkdir(outputDir); end
addpath(fullfile(projectDir,'patient_pinn'),fullfile(projectDir,'synthetic_data'));
S=load(dataFile); ids=unique(pinnResults.predictions.patient_id);
if ~isempty(o.patient_ids); ids=ids(ismember(ids,o.patient_ids)); end
if isempty(ids); error('No requested patients found in PINN results.'); end

profile=[]; bands=[]; coverage=[];
profileFits=cell(numel(ids),1);
for ii=1:numel(ids)
    id=ids(ii); fprintf('\nGamma profile patient %d/%d (ID %d)\n',ii,numel(ids),id);
    center=parameter_value(pinnResults.parameters,id,"Gamma_A","estimate");
    truthGamma=parameter_value(pinnResults.parameters,id,"Gamma_A","truth");
    factors=linspace(1-o.gamma_span,1+o.gamma_span,o.gamma_grid_size);
    grid=unique(max(o.gamma_bounds(1),min(o.gamma_bounds(2),center*factors)));
    fits=cell(numel(grid),1); stat=zeros(numel(grid),1);
    for j=1:numel(grid)
        fprintf('  Gamma %.4f (%d/%d)\n',grid(j),j,numel(grid));
        po=struct('gamma_mode',"population",'gamma_population',grid(j), ...
            'patient_ids',id,'max_patients',1,'max_iterations',o.profile_iterations, ...
            'min_iterations',min(o.profile_min_iterations,o.profile_iterations), ...
            'num_restarts',o.profile_restarts,'refit_iterations',o.profile_refit_iterations, ...
            'use_lbfgs',o.profile_use_lbfgs,'lbfgs_iterations',o.profile_lbfgs_iterations, ...
            'num_collocation',o.num_collocation,'fast_nx',o.fast_nx, ...
            'print_every',o.print_every,'make_plots',false,'save_results',false, ...
            'random_seed',o.random_seed+j);
        fits{j}=run_patient_pinn(dataFile,fullfile(outputDir,'temporary'),po);
        P=fits{j}.predictions; m=P.used_for_fit;
        stat(j)=sum(((P.plateau_pred_ppm(m)-P.plateau_observed_ppm(m))/ ...
            o.acetone_sigma_ppm).^2);
    end
    delta=stat-min(stat); cutoff=2*gammaincinv(o.confidence_level,.5);
    rawAccepted = delta <= cutoff;
    [~,bestj] = min(stat);
    
    % Retain all profile points satisfying the likelihood-ratio threshold.
    accepted = rawAccepted;
    if ~any(accepted); accepted(bestj) = true; end

    lowerTruncated = accepted(1);
    upperTruncated = accepted(end);
    
    profileFits{ii}=fits;
    lo=min(grid(accepted)); hi=max(grid(accepted));
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

sensitivity=dimensionless_sensitivity(pinnResults,S,ids,o);
summary=table(o.confidence_level,height(coverage),mean(coverage.Gamma_covered), ...
    mean(coverage.BHB_holdout_coverage),mean(coverage.slope_holdout_coverage), ...
    mean(coverage.mean_BHB_band_width_mM), ...
    'VariableNames',{'confidence_level','N_patients','Gamma_coverage', ...
    'BHB_holdout_coverage','slope_holdout_coverage','mean_BHB_band_width_mM'});
U=struct('summary',summary,'coverage',coverage,'profile',profile,'bands',bands, ...
    'sensitivity',sensitivity,'options',o,'profile_fits',{profileFits});
writetable(summary,fullfile(outputDir,'uncertainty_summary.csv'));
writetable(coverage,fullfile(outputDir,'uncertainty_coverage.csv'));
writetable(profile,fullfile(outputDir,'gamma_profiles.csv'));
writetable(bands,fullfile(outputDir,'state_uncertainty_bands.csv'));
writetable(sensitivity,fullfile(outputDir,'dimensionless_sensitivity.csv'));
save(fullfile(outputDir,'pinn_uncertainty_results.mat'),'U','-v7.3');
make_plots(U,outputDir); disp(summary); disp(coverage)
end

function T=dimensionless_sensitivity(R,S,ids,o)
T=[]; factors=o.sensitivity_factors(:); base=S.base;
for ii=1:numel(ids)
    id=ids(ii); M=S.meta(S.meta.patient_id==id,:); [~,ix]=sort(M.breath_time_h); M=M(ix,:);
    P=R.predictions(R.predictions.patient_id==id,:); [~,ix]=sort(P.time_h); P=P(ix,:);
    gamma=parameter_value(R.parameters,id,"Gamma_A","estimate");
    a=P.plateau_pred_ppm./max(eps,pulmonary_gain_vector(M,base,o.fast_nx))./ ...
        (base.slow.kappa_QA*base.scales.K_star_mM/gamma/base.scales.lambda_ba*24450);
    baseline=predicted_plateau(M,base,a,gamma,1,1,1,o.fast_nx);
    for name=["Pe","Da_aw","Gamma_A"]
        for j=1:numel(factors)
            fp=1; fd=1; fg=1;
            if name=="Pe"; fp=factors(j); elseif name=="Da_aw"; fd=factors(j); else; fg=factors(j); end
            y=predicted_plateau(M,base,a,gamma,fp,fd,fg,o.fast_nx);
            row=table(id,string(M.scenario(1)),name,factors(j),mean(y), ...
                100*mean((y-baseline)./max(abs(baseline),eps)), ...
                'VariableNames',{'patient_id','scenario','group','factor', ...
                'mean_plateau_ppm','mean_percent_change'});
            T=append(T,row);
        end
    end
end
end

function y=predicted_plateau(M,base,a,gamma,fp,fd,fg,nx)
gain=zeros(height(M),1); p=base; p.fast.nx=nx;
for k=1:height(M)
    p.fast.Pe=M.Pe(k)*fp; p.fast.Da_aw=M.Da_aw(k)*fd;
    z=solve_fast_breath(1,p,M.U_scale(k)); nt=max(2,round(.1*numel(z.time_s)));
    gain(k)=mean(z.c(end,end-nt:end));
end
g=gamma*fg; cstar=base.slow.kappa_QA*base.scales.K_star_mM/g/ ...
    base.scales.lambda_ba*24450; y=gain.*cstar.*a;
end
function g=pulmonary_gain_vector(M,base,nx)
g=zeros(height(M),1); p=base; p.fast.nx=nx;
for k=1:height(M); p.fast.Pe=M.Pe(k); p.fast.Da_aw=M.Da_aw(k);
 z=solve_fast_breath(1,p,M.U_scale(k)); nt=max(2,round(.1*numel(z.time_s)));
 g(k)=mean(z.c(end,end-nt:end)); end
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
d=struct('patient_ids',[],'confidence_level',.95,'gamma_span',.30, ...
 'gamma_grid_size',9,'gamma_bounds',[.3 4.5],'acetone_sigma_ppm',.03, ...
 'profile_iterations',5000,'profile_min_iterations',2500,'profile_restarts',1, ...
 'profile_refit_iterations',500,'profile_use_lbfgs',true, ...
 'profile_lbfgs_iterations',100,'num_collocation',121,'fast_nx',81, ...
 'print_every',1000,'random_seed',73001,'sensitivity_factors',[.5 .75 1 1.25 1.5], ...
 'monitor_bhb_mM',1.5,'urgent_bhb_mM',3.0,'rising_mM_h',.15, ...
 'rapid_rise_mM_h',.30);
f=fieldnames(d); for k=1:numel(f); if ~isfield(o,f{k}); o.(f{k})=d.(f{k}); end; end
end
