function report = run_recoverability_analysis(cohort,outputDir)
%RUN_RECOVERABILITY_ANALYSIS Test rising/falling identifiability.
% Finds developing/resolving breaths with nearly equal instantaneous breath
% acetone, then measures whether recent breath history separates them.
if nargin<2 || isempty(outputDir)
 outputDir=fullfile(fileparts(mfilename('fullpath')),'outputs','cohort');
end
if ischar(cohort) || isstring(cohort)
 S=load(cohort,'meta'); meta=S.meta;
elseif isstruct(cohort) && isfield(cohort,'metadata'); meta=cohort.metadata;
else; meta=cohort; end

use=meta.scenario~="baseline";
M=meta(use,:); M.direction=double(M.dBHB_dt_mM_h>0);
dev=find(M.direction==1); res=find(M.direction==0);
used=false(numel(res),1); pairs=zeros(0,2); mismatch=[];
for i=1:numel(dev)
 [e,j]=min(abs(M.C_exh_plateau_ppm(dev(i))-M.C_exh_plateau_ppm(res)));
 if ~used(j)
   used(j)=true; pairs(end+1,:)=[dev(i) res(j)]; mismatch(end+1,1)=e; %#ok<AGROW>
 end
end
P1=M(pairs(:,1),:); P0=M(pairs(:,2),:);
pairTable=table(P1.patient_id,P1.breath_time_h,P1.BHB_mM,P1.dBHB_dt_mM_h, ...
 P0.patient_id,P0.breath_time_h,P0.BHB_mM,P0.dBHB_dt_mM_h,mismatch, ...
 'VariableNames',{'rising_patient','rising_time_h','rising_BHB_mM','rising_slope', ...
 'falling_patient','falling_time_h','falling_BHB_mM','falling_slope','plateau_mismatch_ppm'});
writetable(pairTable,fullfile(outputDir,'matched_rising_falling_pairs.csv'));

% A simple, toolbox-free direction score: change from the patient's prior breath.
M.delta_plateau=nan(height(M),1); M.delta_auc=nan(height(M),1);
for id=unique(M.patient_id).'
 idx=find(M.patient_id==id); [~,ord]=sort(M.breath_time_h(idx)); idx=idx(ord);
 M.delta_plateau(idx(2:end))=diff(M.C_exh_plateau_ppm(idx));
 M.delta_auc(idx(2:end))=diff(M.C_exh_auc_ppm_s(idx));
end
valid=isfinite(M.delta_plateau);
predLong=M.delta_plateau(valid)>0;
accuracyLong=mean(predLong==M.direction(valid));
% Instantaneous plateau threshold selected on training data, evaluated on test.
tr=string(M.split)=="train"; te=string(M.split)=="test";
thresholds=linspace(min(M.C_exh_plateau_ppm(tr)),max(M.C_exh_plateau_ppm(tr)),200);
accHigh=zeros(size(thresholds)); accLow=accHigh;
for k=1:numel(thresholds)
 accHigh(k)=mean((M.C_exh_plateau_ppm(tr)>thresholds(k))==M.direction(tr));
 accLow(k)=mean((M.C_exh_plateau_ppm(tr)<thresholds(k))==M.direction(tr));
end
[bestHigh,kHigh]=max(accHigh); [bestLow,kLow]=max(accLow);
if bestHigh>=bestLow
 threshold=thresholds(kHigh); instantSign=1;
else
 threshold=thresholds(kLow); instantSign=-1;
end
accuracyInstant=mean(((instantSign*M.C_exh_plateau_ppm(te))> ...
    instantSign*threshold)==M.direction(te));

report=table(height(pairTable),median(mismatch),accuracyInstant,accuracyLong,threshold, ...
 'VariableNames',{'matched_pairs','median_plateau_mismatch_ppm', ...
 'instantaneous_test_accuracy','longitudinal_accuracy','instantaneous_threshold_ppm'});
writetable(report,fullfile(outputDir,'recoverability_summary.csv'));
save(fullfile(outputDir,'recoverability_results.mat'),'report','pairTable','M');

fig=figure('Color','w','Position',[100 100 1000 420]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile; groupedScatter(M.C_exh_plateau_ppm,M.BHB_mM,M.scenario); grid on
xlabel('Exhaled acetone plateau (ppm)'); ylabel('BHB (mmol/L)');
title('Same breath acetone can imply different BHB states');
nexttile; groupedScatter(M.delta_plateau,M.dBHB_dt_mM_h,M.scenario); grid on
xline(0,'k--'); yline(0,'k--'); xlabel('\Delta breath acetone (ppm)');
ylabel('dBHB/dt (mmol L^{-1} h^{-1})'); title('History recovers direction');
exportgraphics(fig,fullfile(outputDir,'recoverability.png'),'Resolution',220);
disp(report)
end

function groupedScatter(x,y,group)
% No Statistics and Machine Learning Toolbox required.
hold on; names=unique(string(group)); colors=lines(numel(names));
for k=1:numel(names)
    use=string(group)==names(k) & isfinite(x) & isfinite(y);
    scatter(x(use),y(use),22,colors(k,:),'filled','MarkerFaceAlpha',0.65, ...
        'DisplayName',names(k));
end
legend('Location','best');
end
