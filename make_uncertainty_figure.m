%% MAKE_UNCERTAINTY_FIGURE
% Top row: profile likelihoods for representative stable and developing cases.
% Bottom row: BHB point estimate, profile envelope, truth, and alert thresholds.

clear; close all; clc;

% -------------------------------------------------------------------------
% Locate the saved uncertainty-analysis output.
% -------------------------------------------------------------------------
root=fileparts(mfilename('fullpath'));
resultFile=fullfile(root,'pinn_uncertainty','outputs','full','pinn_uncertainty_results.mat');
assert(isfile(resultFile),'Missing profile results.');
fprintf('Using uncertainty results: %s\n',resultFile);

S = load(resultFile,'U');
U = S.U;

P = U.profile;
B = U.bands;
C = U.coverage;

% -------------------------------------------------------------------------
% Select representative patients.
% Choose the stable patient with the broadest accepted Gamma_A interval
% and the developing patient with the narrowest accepted interval.
% -------------------------------------------------------------------------
C.gamma_width = C.Gamma_upper - C.Gamma_lower;

stableRows = find(string(C.scenario) == "baseline");
developingRows = find(string(C.scenario) == "developing");

if isempty(stableRows) || isempty(developingRows)
    error('The uncertainty results must contain both baseline and developing patients.');
end

[~,stableIndex] = max(C.gamma_width(stableRows));
[~,developingIndex] = min(C.gamma_width(developingRows));

stableID = C.patient_id(stableRows(stableIndex));
developingID = C.patient_id(developingRows(developingIndex));

ids = [stableID, developingID];
labels = ["Stable", "Developing"];

fprintf('Selected stable patient %d and developing patient %d.\n', ...
    stableID,developingID);

% -------------------------------------------------------------------------
% Visual settings.
% -------------------------------------------------------------------------
profileColor  = [0.00 0.45 0.70];
bandColor     = [0.72 0.85 0.95];
truthColor    = [0.10 0.10 0.10];
pointColor    = [0.00 0.45 0.70];
monitorColor  = [0.90 0.50 0.10];
urgentColor   = [0.75 0.20 0.15];

fig = figure( ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[1 1 7.2 5.2], ...
    'Renderer','painters');

set(fig, ...
    'DefaultAxesFontName','Times New Roman', ...
    'DefaultTextFontName','Times New Roman');

tiledlayout(2,2,'TileSpacing','loose','Padding','compact');

% -------------------------------------------------------------------------
% Top row: Gamma_A profile likelihoods.
% -------------------------------------------------------------------------
for k = 1:2
    Q = P(P.patient_id == ids(k),:);
    Q = sortrows(Q,'Gamma_A');

    nexttile;
    hold on;

    hProfile = plot(Q.Gamma_A,Q.delta_chi_square, ...
        '-o', ...
        'Color',profileColor, ...
        'MarkerFaceColor','w', ...
        'MarkerSize',4.5, ...
        'LineWidth',1.4);

    hAccepted = scatter(Q.Gamma_A(Q.accepted), ...
        Q.delta_chi_square(Q.accepted), ...
        26,profileColor,'filled');

    hCutoff = yline(Q.cutoff(1),'k--','LineWidth',0.9);

    truthGamma = C.Gamma_truth(C.patient_id == ids(k));
    hTruth = xline(truthGamma,'-', ...
        'Color',truthColor, ...
        'LineWidth',1.0);

    grid on;
    box on;

    xlabel('\Gamma_A','Interpreter','tex');
    ylabel('\Delta\chi^2','Interpreter','tex');

    title(sprintf('%s patient %d: parameter profile', ...
        labels(k),ids(k)), ...
        'FontName','Times New Roman', ...
        'FontSize',10, ...
        'FontWeight','bold');

    set(gca, ...
        'FontName','Times New Roman', ...
        'FontSize',8.5, ...
        'LineWidth',0.7);

    % A compact legend appears only once.
    if k == 1
        legend([hProfile,hAccepted,hTruth,hCutoff], ...
            {'Profile fit','Accepted refit','True \Gamma_A','Retention cutoff'}, ...
            'Location','northeast', ...
            'FontName','Times New Roman', ...
            'FontSize',6);
    end
end

% -------------------------------------------------------------------------
% Bottom row: BHB profile envelopes and research alert thresholds.
% -------------------------------------------------------------------------
for k = 1:2
    Q = B(B.patient_id == ids(k),:);
    Q = sortrows(Q,'time_h');

    nexttile;
    hold on;

    hBand = fill([Q.time_h; flipud(Q.time_h)], ...
        [Q.BHB_lower_mM; flipud(Q.BHB_upper_mM)], ...
        bandColor, ...
        'EdgeColor','none', ...
        'FaceAlpha',0.85);

    hTruth = plot(Q.time_h,Q.BHB_truth_mM, ...
        '-','Color',truthColor,'LineWidth',1.6);

    hPoint = plot(Q.time_h,Q.BHB_point_mM, ...
        '-','Color',pointColor,'LineWidth',1.8);

    % No automatic labels: they crowded the original figure.
    yline(1.5,'--', ...
        'Color',monitorColor, ...
        'LineWidth',1.0);

    yline(3.0,'--', ...
        'Color',urgentColor, ...
        'LineWidth',1.0);

    % Keep both BHB alert thresholds visible.
    ymax = max([Q.BHB_upper_mM; Q.BHB_truth_mM; 3.0]);
    if k == 1
        ylim([min(Q.BHB_truth_mM)*0.9, max(Q.BHB_upper_mM)*1.05]);
    else
        ylim([0,1.06*ymax]);
    end

    grid on;
    box on;

    xlabel('Time (h)');
    ylabel('BHB (mmol/L)');

    title(sprintf('%s patient %d: BHB profile envelope', ...
        labels(k),ids(k)), ...
        'FontName','Times New Roman', ...
        'FontSize',10, ...
        'FontWeight','bold');

    set(gca, ...
        'FontName','Times New Roman', ...
        'FontSize',8.5, ...
        'LineWidth',0.7);

    % A compact legend appears only once.
    if k == 1
        legend([hBand,hTruth,hPoint], ...
            {'Profile sensitivity envelope','Simulated BHB','Minimum-misfit profile refit'}, ...
            'Location','northeast', ...
            'FontName','Times New Roman', ...
            'FontSize',6);
    end
end

% -------------------------------------------------------------------------
% Export.
% -------------------------------------------------------------------------
outputDir = fileparts(resultFile);

exportgraphics(fig,fullfile(outputDir,'uncertainty_figure.pdf'), ...
    'ContentType','vector');

exportgraphics(fig,fullfile(outputDir,'uncertainty_figure.png'), ...
    'Resolution',600);

fprintf('\nSaved:\n%s\n%s\n', ...
    fullfile(outputDir,'uncertainty_figure.pdf'), ...
    fullfile(outputDir,'uncertainty_figure.png'));