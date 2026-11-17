%% MAKE_HISTORY_FIGURE
% Creates the two-panel figure: Why breath history matters.
% Run this script from the top-level folder of your project.

clear; close all; clc;

root=fileparts(mfilename('fullpath'));
metadataFile=fullfile(root,'synthetic_data','outputs','cohort','cohort_metadata.csv');
assert(isfile(metadataFile),'Missing cohort metadata.');
fprintf('Using cohort metadata: %s\n', metadataFile);

M = readtable(metadataFile);

% Keep only developing and resolving virtual patients.
isDev = string(M.scenario) == "developing";
isRes = string(M.scenario) == "resolving";
use = isDev | isRes;
M = M(use,:);

% Compute change in plateau acetone from each patient's previous breath.
M.delta_plateau_ppm = nan(height(M),1);

patientIDs = unique(M.patient_id);
for k = 1:numel(patientIDs)
    rows = find(M.patient_id == patientIDs(k));
    [~,order] = sort(M.breath_time_h(rows));
    rows = rows(order);

    M.delta_plateau_ppm(rows(2:end)) = ...
        diff(M.C_exh_plateau_ppm(rows));
end

isDev = string(M.scenario) == "developing";
isRes = string(M.scenario) == "resolving";
hasHistory = isfinite(M.delta_plateau_ppm);

% Print-friendly colors.
blue   = [0.00 0.45 0.70];
orange = [0.85 0.33 0.10];

fig = figure( ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[1 1 7.2 3.15], ...
    'Renderer','painters');

tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Panel A: Acetone alone is ambiguous.
nexttile;
hold on;

scatter(M.C_exh_plateau_ppm(isDev), M.BHB_mM(isDev), ...
    25, blue, 'filled', ...
    'MarkerFaceAlpha',0.55, ...
    'DisplayName','Developing ketosis');

scatter(M.C_exh_plateau_ppm(isRes), M.BHB_mM(isRes), ...
    25, orange, 'filled', ...
    'MarkerFaceAlpha',0.55, ...
    'DisplayName','Resolving ketosis');

grid on;
box on;
xlabel('Breath-acetone plateau (ppm)');
ylabel('Blood BHB (mmol/L)');
t1 = title('Instantaneous acetone does not uniquely determine BHB');
t1.Units = 'normalized';
t1.Position(2) = 1.05;
legend('Location','northwest','FontSize',8);
set(gca,'FontName','Times New Roman','FontSize',9.5,'LineWidth',0.75);

% Panel B: Change from prior breath reveals direction.
nexttile;
hold on;

scatter(M.delta_plateau_ppm(isDev & hasHistory), ...
        M.dBHB_dt_mM_h(isDev & hasHistory), ...
    25, blue, 'filled', ...
    'MarkerFaceAlpha',0.55);

scatter(M.delta_plateau_ppm(isRes & hasHistory), ...
        M.dBHB_dt_mM_h(isRes & hasHistory), ...
    25, orange, 'filled', ...
    'MarkerFaceAlpha',0.55);

xline(0,'k--','LineWidth',0.8);
yline(0,'k--','LineWidth',0.8);

grid on;
box on;
xlabel('Change from preceding breath (ppm)');
ylabel('BHB slope (mmol/L/h)');
t2 = title('Breath history helps recover BHB direction');
t2.Units = 'normalized';
t2.Position(2) = 1.05;
set(gca,'FontName','Times New Roman','FontSize',9.5,'LineWidth',0.75);

% Save beside the metadata file.
outputDir = fileparts(metadataFile);

exportgraphics(fig,fullfile(outputDir,'history_figure.pdf'), ...
    'ContentType','vector');

exportgraphics(fig,fullfile(outputDir,'history_figure.png'), ...
    'Resolution',600);

fprintf('\nSaved figure files:\n%s\n%s\n', ...
    fullfile(outputDir,'history_figure.pdf'), ...
    fullfile(outputDir,'history_figure.png'));