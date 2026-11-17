function cohort = generate_synthetic_cohort(outputDir,opts)
%GENERATE_SYNTHETIC_COHORT Create noisy breath spectra and hidden truth.
% One observation is one breath. Raw input arrays have dimensions
% breath-time x wavelength x observation. Fixed random seeds make the
% dataset exactly reproducible.

if nargin<1 || isempty(outputDir)
    outputDir=fullfile(fileparts(mfilename('fullpath')),'outputs','cohort');
end
if nargin<2; opts=struct; end
opts=defaults(opts);
if ~exist(outputDir,'dir'); mkdir(outputDir); end
rng(opts.seed,'twister');

base=ketone_parameters();
scenarios={'baseline','developing','resolving'};
nObs=opts.patients_per_scenario*numel(scenarios)*numel(opts.breath_hours);
nT=opts.breath_samples; nW=opts.wavelength_samples;
spectra=zeros(nT,nW,nObs,'single'); cleanSpectra=spectra;
time_s=zeros(nT,nObs,'single');
meta=table('Size',[nObs 21], ...
 'VariableTypes',{'double','string','double','double','double','double', ...
 'double','double','double','double','double','double','double','double', ...
 'double','double','double','double','double','double','double'}, ...
 'VariableNames',{'patient_id','scenario','breath_time_h','BHB_mM', ...
 'dBHB_dt_mM_h','AcAc_mM','blood_acetone_uM','a_dimensionless','r0', ...
 'delta_r','tau_c','k_r','Pe','Da_aw','Gamma_A','A_star_mM','C_star_ppm','U_scale', ...
 'C_exh_plateau_ppm','C_exh_auc_ppm_s','noise_sd_abs'});

obs=0; patient=0; wn=[];
for is=1:numel(scenarios)
    scenario=scenarios{is};
    for ip=1:opts.patients_per_scenario
        patient=patient+1; p=sample_patient(base,scenario,opts);
        slow=solve_slow_metabolism(p,scenario);
        for ib=1:numel(opts.breath_hours)
            obs=obs+1; h=opts.breath_hours(ib);
            b=interp1(slow.time_h,slow.b,h,'pchip');
            q=interp1(slow.time_h,slow.q,h,'pchip');
            a=interp1(slow.time_h,slow.a,h,'pchip');
            r=interp1(slow.time_h,slow.r,h,'pchip');
            db=p.slow.rho_B*r-(1+p.slow.kappa_BQ)*b+p.slow.kappa_QB*q;
            U=max(opts.U_bounds(1),min(opts.U_bounds(2),1+opts.U_cv*randn));
            breath=solve_fast_breath(a,p,U);
            tq=linspace(0,p.fast.breath_duration_s,nT).';
            Cq=interp1(breath.time_s,breath.C_exh_ppm,tq,'pchip');
            [wnFull,Afull]=synthetic_acetone_spectrum(Cq,p);
            idx=round(linspace(1,numel(wnFull),nW)); wn=wnFull(idx);
            Aclean=Afull(:,idx);
            sigma=max(opts.noise_floor_abs,opts.noise_relative*max(Aclean(:)));
            Anoise=Aclean+sigma*randn(size(Aclean));
            cleanSpectra(:,:,obs)=single(Aclean);
            spectra(:,:,obs)=single(Anoise);
            time_s(:,obs)=single(tq);
            meta(obs,:)={patient,string(scenario),h,p.scales.K_star_mM*b, ...
                p.scales.K_star_mM/p.scales.T_m_h*db, ...
                p.scales.K_star_mM*q,1000*p.scales.A_star_mM*a,a, ...
                p.slow.r0,p.slow.delta_r,p.slow.tau_c,p.slow.k_r, ...
                p.fast.Pe,p.fast.Da_aw,p.slow.Gamma_A, ...
                p.scales.A_star_mM,p.scales.C_star_ppm,U, ...
                mean(Cq(end-max(2,round(.1*nT)):end)),trapz(tq,Cq),sigma};
        end
    end
end

% Stratified patient-level split prevents leakage and retains each scenario.
split=repmat("test",height(meta),1);
ids=unique(meta.patient_id);
for is=1:numel(scenarios)
    scenarioIds=unique(meta.patient_id(meta.scenario==string(scenarios{is})));
    scenarioIds=scenarioIds(randperm(numel(scenarioIds)));
    nTrain=max(1,round(opts.train_fraction*numel(scenarioIds)));
    nVal=max(0,round(opts.validation_fraction*numel(scenarioIds)));
    nVal=min(nVal,numel(scenarioIds)-nTrain);
    split(ismember(meta.patient_id,scenarioIds(1:nTrain)))="train";
    split(ismember(meta.patient_id,scenarioIds(nTrain+1:nTrain+nVal)))="validation";
end
meta.split=categorical(split,{'train','validation','test'});
writetable(meta,fullfile(outputDir,'cohort_metadata.csv'));
cohort_version=2;
save(fullfile(outputDir,'synthetic_cohort.mat'),'spectra','cleanSpectra', ...
    'time_s','wn','meta','opts','base','cohort_version','-v7.3');
cohort=struct('spectra',spectra,'cleanSpectra',cleanSpectra, ...
    'time_s',time_s,'wavenumber_cm1',wn,'metadata',meta,'options',opts, ...
    'version',cohort_version);
fprintf('Saved %d breaths from %d virtual patients to %s\n', ...
    height(meta),numel(ids),outputDir);
end

function p=sample_patient(p,scenario,o)
p.slow.r0=uniform(o.r0_range);
p.slow.delta_r=uniform(o.delta_r_range);
p.slow.tau_c=uniform(o.tau_c_range);
p.slow.k_r=uniform(o.k_r_range);
p.fast.Pe=loguniform(o.Pe_range);
p.fast.Da_aw=loguniform(o.Da_aw_range);
p.slow.Gamma_A=loguniform(o.Gamma_A_range);
% A_* and C_* are patient-specific because Gamma_A=k_A,eff/k_B.
p.scales.A_star_mM=p.slow.kappa_QA*p.scales.K_star_mM/p.slow.Gamma_A;
p.scales.C_star_ppm=p.scales.A_star_mM/p.scales.lambda_ba*24450;
if strcmp(scenario,'baseline'); p.slow.delta_r=0; end
end
function x=uniform(r); x=r(1)+(r(2)-r(1))*rand; end
function x=loguniform(r); x=exp(log(r(1))+(log(r(2))-log(r(1)))*rand); end
function o=defaults(o)
d=struct('seed',271828,'patients_per_scenario',12, ...
 'breath_hours',0.5:0.75:8,'breath_samples',101,'wavelength_samples',61, ...
 'r0_range',[0.15 0.30],'delta_r_range',[1.5 2.7], ...
 'tau_c_range',[2.3 3.7],'k_r_range',[1.5 3.0], ...
 'Pe_range',[12 60],'Da_aw_range',[0.35 1.5], ...
 'Gamma_A_range',[0.7 2.8],'U_cv',0.10,'U_bounds',[0.7 1.3], ...
 'noise_relative',0.02,'noise_floor_abs',2e-6, ...
 'train_fraction',0.65,'validation_fraction',0.15);
f=fieldnames(d); for k=1:numel(f); if ~isfield(o,f{k}); o.(f{k})=d.(f{k}); end; end
end
