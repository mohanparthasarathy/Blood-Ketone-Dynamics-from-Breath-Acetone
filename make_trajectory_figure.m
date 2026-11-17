%% MAKE_TRAJECTORY_FIGURE
% Six held-out virtual patients: two stable, two developing, two resolving.
% Uses the practical prior-mode PINN result when it finds one.

clear; close all; clc;

root=fileparts(mfilename('fullpath'));
resultFile=fullfile(root,'patient_pinn','outputs','prior10_v2','patient_pinn_results.mat');
assert(isfile(resultFile),'Missing informed-prior results.');
fprintf('Using PINN results: %s\n',resultFile);

S = load(resultFile,'results');

if ~isfield(S,'results') || ~isfield(S.results,'predictions')
    error('The selected MAT-file does not contain results.predictions.');
end

P = S.results.predictions;

required = {'patient_id','time_h','scenario','used_for_fit', ...
    'temporal_holdout','BHB_true_mM','BHB_pred_mM'};

if ~all(ismember(required,P.Properties.VariableNames))
    error('The selected prediction table does not contain the expected variables.');
end

% Order panels as stable, developing, resolving.
scenarioOrder = ["baseline","developing","resolving"];
displayNames = ["Stable","Developing","Resolving"];

patientIDs = [];
patientScenarios = strings(0);

for k = 1:numel(scenarioOrder)
    ids = unique(P.patient_id(string(P.scenario)==scenarioOrder(k)));
    ids = ids(:).';

    if numel(ids) < 2
        error('Expected at least two held-out patients in scenario "%s".', ...
            scenarioOrder(k));
    end

    patientIDs = [patientIDs, ids(1:2)]; %#ok<AGROW>
    patientScenarios = [patientScenarios, ...
        repmat(displayNames(k),1,2)]; %#ok<AGROW>
end

% Use one common vertical scale so the six metabolic states are comparable.
yMax = 1.06 * max([P.BHB_true_mM; P.BHB_pred_mM]);
yMin = 0;
yRange = yMax - yMin;
markerY = yMin + 0.045*yRange;

truthColor = [0.10 0.10 0.10];
pinnColor  = [0.00 0.45 0.70];
fitColor   = [0.25 0.25 0.25];
holdColor  = [0.80 0.20 0.15];

fig = figure( ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[1 1 7.2 4.9], ...
    'Renderer','painters');

tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

for k = 1:numel(patientIDs)
    Q = P(P.patient_id==patientIDs(k),:);
    Q = sortrows(Q,'time_h');

    nexttile;
    hold on;

    % True and reconstructed BHB trajectories.
    hTruth = plot(Q.time_h,Q.BHB_true_mM, ...
        '-','Color',truthColor,'LineWidth',1.6);

    hPINN = plot(Q.time_h,Q.BHB_pred_mM, ...
        '-','Color',pinnColor,'LineWidth',1.8);

    % Small markers at the bottom identify measurement use without
    % implying that BHB itself was measured at those times.
    fitTimes = Q.time_h(Q.used_for_fit);
    holdTimes = Q.time_h(Q.temporal_holdout);

    hFit = scatter(fitTimes,markerY*ones(size(fitTimes)), ...
        18,fitColor,'filled','Marker','o');

    hHold = scatter(holdTimes,markerY*ones(size(holdTimes)), ...
        34,holdColor,'filled','Marker','v');

    grid on;
    box on;
    xlim([min(Q.time_h)-0.15, max(Q.time_h)+0.15]);
    ylim([yMin,yMax]);

    % Label only the outer axes to avoid crowding.
    if k > 3
        xlabel('Time (h)');
    else
        xlabel('');
    end
    
    if mod(k-1,3) == 0
        ylabel('BHB (mmol/L)');
    else
        ylabel('');
    end

    panelTitle = sprintf('%s patient %d', ...
        patientScenarios(k),patientIDs(k));
    
    title(panelTitle,'FontWeight','bold');

    set(gca,'FontName','Times New Roman', ...
        'FontSize',8.5,'LineWidth',0.7);
    if k == 1
        legend([hTruth,hPINN,hFit,hHold], ...
            {'Simulated BHB','PINN estimate', ...
             'Fit breath','Withheld breath'}, ...
            'Location','northwest', ...
            'FontName','Times New Roman', ...
            'FontSize',6.5);
    end
end


% Save beside the results file.
outputDir = fileparts(resultFile);

exportgraphics(fig,fullfile(outputDir,'trajectory_figure.pdf'), ...
    'ContentType','vector');

exportgraphics(fig,fullfile(outputDir,'trajectory_figure.png'), ...
    'Resolution',600);

fprintf('\nSaved:\n%s\n%s\n', ...
    fullfile(outputDir,'trajectory_figure.pdf'), ...
    fullfile(outputDir,'trajectory_figure.png'));