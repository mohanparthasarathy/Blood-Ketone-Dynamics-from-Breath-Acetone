function stats=run_fisher_validation_statistics(resultFile,outputDir,userOpts)
%RUN_FISHER_VALIDATION_STATISTICS Patient-clustered paired analysis.
% Bootstrap resamples patients, preserving all holdout breaths within patient.

here=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(resultFile)
    resultFile=fullfile(here,'outputs','fisher_validation_full', ...
        'fisher_pinn_validation.mat');
end
if nargin<2 || isempty(outputDir)
    outputDir=fullfile(here,'outputs','fisher_validation_statistics');
end
if nargin<3; userOpts=struct; end
o=defaults(userOpts);
if ~isfile(resultFile); error('Results file not found: %s',resultFile); end
if ~exist(outputDir,'dir'); mkdir(outputDir); end

S=load(resultFile);
if isfield(S,'validation'); V=S.validation;
else; error('The MAT file does not contain variable validation.'); end
R=V.results; labels=string(V.summary.design); nD=numel(R);
ids=unique(R{1}.predictions.patient_id).'; nPatients=numel(ids);
metrics=["BHB_RMSE_mM","slope_RMSE_mM_h","acetone_RMSE_ppm","Gamma_A_MAE"];
nM=numel(metrics); point=zeros(nD,nM); patientValues=zeros(nPatients,nD,nM);

for d=1:nD
    point(d,:)=calculate_metrics(R{d},ids);
    for i=1:nPatients
        patientValues(i,d,:)=calculate_metrics(R{d},ids(i));
    end
end

% Cluster bootstrap: resample complete patients with replacement.
rng(o.random_seed,'twister'); boot=zeros(o.bootstrap_replicates,nD,nM);
for b=1:o.bootstrap_replicates
    sampled=ids(randi(nPatients,1,nPatients));
    for d=1:nD; boot(b,d,:)=calculate_metrics(R{d},sampled); end
end

alpha=(1-o.confidence_level)/2;
absolute=cell(0,1+3*nM);
for d=1:nD
    row={labels(d)};
    for m=1:nM
        ci=empirical_ci(squeeze(boot(:,d,m)),alpha);
        row=[row,{point(d,m),ci(1),ci(2)}]; %#ok<AGROW>
    end
    absolute=[absolute;row]; %#ok<AGROW>
end
absolute=cell2table(absolute,'VariableNames',absolute_names(metrics));
% for j=2:width(absolute); absolute{:,j}=cell2mat(absolute{:,j}); end

% Fisher-minus-comparator differences. Negative means Fisher has less error.
fisher=find(labels=="fisher_D8",1); competitors=find((1:nD)~=fisher);
paired=[];
for c=competitors
    for m=1:nM
        delta=point(fisher,m)-point(c,m);
        bdelta=squeeze(boot(:,fisher,m)-boot(:,c,m)); ci=empirical_ci(bdelta,alpha);
        fv=squeeze(patientValues(:,fisher,m)); cv=squeeze(patientValues(:,c,m));
        wins=sum(fv<cv); ties=sum(abs(fv-cv)<=o.tie_tolerance);
        row=table(labels(c),metrics(m),point(fisher,m),point(c,m),delta, ...
            100*delta/max(eps,point(c,m)),ci(1),ci(2),wins,ties,nPatients, ...
            'VariableNames',{'comparator','metric','Fisher_value', ...
            'comparator_value','Fisher_minus_comparator','percent_difference', ...
            'CI_lower','CI_upper','Fisher_patient_wins','ties','N_patients'});
        paired=append_table(paired,row);
    end
end

patientTable=[];
for i=1:nPatients
    for d=1:nD
        row=table(ids(i),labels(d),patientValues(i,d,1),patientValues(i,d,2), ...
            patientValues(i,d,3),patientValues(i,d,4), ...
            'VariableNames',{'patient_id','design','BHB_RMSE_mM', ...
            'slope_RMSE_mM_h','acetone_RMSE_ppm','Gamma_A_abs_error'});
        patientTable=append_table(patientTable,row);
    end
end

counts=V.summary.measurement_count;
retention=table(labels,counts,100*(1-counts/counts(1)), ...
    point(:,1)/point(1,1),point(:,2)/point(1,2),point(:,3)/point(1,3), ...
    'VariableNames',{'design','measurement_count','measurement_reduction_percent', ...
    'BHB_RMSE_ratio_to_full','slope_RMSE_ratio_to_full', ...
    'acetone_RMSE_ratio_to_full'});

stats=struct('confidence_level',o.confidence_level,'N_patients',nPatients, ...
    'absolute_metrics',absolute,'paired_differences',paired, ...
    'patient_metrics',patientTable,'retention',retention, ...
    'bootstrap_metrics',boot,'options',o);
writetable(absolute,fullfile(outputDir,'bootstrap_absolute_metrics.csv'));
writetable(paired,fullfile(outputDir,'paired_fisher_comparisons.csv'));
writetable(patientTable,fullfile(outputDir,'patient_level_metrics.csv'));
writetable(retention,fullfile(outputDir,'performance_retention.csv'));
save(fullfile(outputDir,'fisher_validation_statistics.mat'),'stats','-v7.3');
make_plot(labels,point,boot,alpha,outputDir);
disp(retention); disp(paired)
end

function v=calculate_metrics(R,sampledIds)
eb=[]; es=[]; ea=[]; eg=[];
for id=sampledIds
    P=R.predictions(R.predictions.patient_id==id & ...
        R.predictions.temporal_holdout,:);
    eb=[eb;P.BHB_pred_mM-P.BHB_true_mM]; %#ok<AGROW>
    es=[es;P.slope_pred_mM_h-P.slope_true_mM_h]; %#ok<AGROW>
    ea=[ea;P.plateau_pred_ppm-P.plateau_observed_ppm]; %#ok<AGROW>
    Q=R.parameters(R.parameters.patient_id==id & R.parameters.parameter=="Gamma_A",:);
    eg=[eg;Q.error]; %#ok<AGROW>
end
v=[sqrt(mean(eb.^2)),sqrt(mean(es.^2)),sqrt(mean(ea.^2)),mean(abs(eg))];
end

function names=absolute_names(metrics)
names={'design'};
for m=1:numel(metrics)
    root=char(metrics(m)); names=[names,{root,[root '_CI_lower'],[root '_CI_upper']}]; %#ok<AGROW>
end
end

function ci=empirical_ci(x,alpha)
x=sort(x(isfinite(x))); n=numel(x);
ci=[interp_quantile(x,alpha,n),interp_quantile(x,1-alpha,n)];
end
function q=interp_quantile(x,p,n)
z=1+(n-1)*p; lo=floor(z); hi=ceil(z);
if lo==hi; q=x(lo); else; q=x(lo)+(z-lo)*(x(hi)-x(lo)); end
end

function make_plot(labels,point,boot,alpha,out)
f=figure('Color','w','Position',[100 100 1250 410]);
titles={'BHB RMSE','Slope RMSE','Acetone RMSE'};
units={'mmol/L','mmol L^{-1} h^{-1}','ppm'};
tiledlayout(1,3,'TileSpacing','compact');
for m=1:3
    lo=zeros(numel(labels),1); hi=lo;
    for d=1:numel(labels)
        ci=empirical_ci(squeeze(boot(:,d,m)),alpha); lo(d)=ci(1); hi(d)=ci(2);
    end
    cats=categorical(labels,labels,'Ordinal',true);
    nexttile; bar(cats,point(:,m)); hold on
    errorbar(1:numel(labels),point(:,m),point(:,m)-lo,hi-point(:,m), ...
        'k.','LineWidth',1.2); grid on; title(titles{m}); ylabel(units{m});
    xtickangle(25)
end
exportgraphics(f,fullfile(out,'bootstrap_design_comparison.png'),'Resolution',220);
close(f)
end

function T=append_table(T,A)
if isempty(T); T=A; else; T=[T;A]; end
end

function o=defaults(o)
d=struct('bootstrap_replicates',10000,'confidence_level',.95, ...
    'random_seed',9132026,'tie_tolerance',1e-12);
f=fieldnames(d);
for k=1:numel(f); if ~isfield(o,f{k}); o.(f{k})=d.(f{k}); end; end
end
