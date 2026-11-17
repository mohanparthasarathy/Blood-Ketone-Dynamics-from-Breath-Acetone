function results = run_ketone_forward_model(outputDir)
%RUN_KETONE_FORWARD_MODEL Run the two-timescale ketone/breath simulator.
%   RESULTS = RUN_KETONE_FORWARD_MODEL() solves the slow metabolic ODE,
%   freezes blood acetone at selected metabolic snapshots, and solves a fast
%   method-of-lines airway model for each breath.

if nargin < 1
    outputDir = fullfile(fileparts(mfilename('fullpath')), 'outputs');
end
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

p = ketone_parameters();
slow = solve_slow_metabolism(p, 'developing');

breathHours = [0.5 3.0 5.0 8.0];
aSamples = interp1(slow.time_h, slow.a, breathHours, 'pchip');
breaths = cell(size(breathHours));
for k = 1:numel(breathHours)
    breaths{k} = solve_fast_breath(aSamples(k), p);
end

selected = breaths{3}; % five-hour metabolic snapshot
[wn, absorbance] = synthetic_acetone_spectrum(selected.C_exh_ppm, p);

fig = figure('Color','w','Position',[80 80 1200 780]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
yyaxis left
plot(slow.time_h,slow.B_mM,'LineWidth',1.8); hold on
plot(slow.time_h,slow.Q_mM,'LineWidth',1.8);
ylabel('Concentration (mmol/L)');
yyaxis right
plot(slow.time_h,slow.r,'--','Color',[0.45 0.45 0.45],'LineWidth',1.5);
ylabel('Dimensionless ketogenesis drive');
xlabel('Metabolic time (h)'); title('Slow ketone dynamics');
legend('BHB','AcAc','Ketogenesis drive','Location','northwest'); grid on

nexttile;
plot(slow.time_h,1000*slow.A_blood_mM,'LineWidth',1.8); hold on
scatter(breathHours,1000*p.scales.A_star_mM*aSamples,32,'k','filled');
xlabel('Metabolic time (h)'); ylabel('Blood acetone (\mumol/L)');
title('Delayed acetone response'); legend('Blood acetone','Breath snapshots'); grid on

nexttile;
profileTimes = [0.25 0.75 1.5 3.0 5.0];
hold on
for k = 1:numel(profileTimes)
    [~,j] = min(abs(selected.time_s-profileTimes(k)));
    plot(selected.x_m,selected.C_ppm(:,j),'LineWidth',1.5, ...
        'DisplayName',sprintf('%.2g s',selected.time_s(j)));
end
xlabel('Effective airway position (m)'); ylabel('Acetone (ppm)');
title('Concentration-front propagation'); legend('Location','best'); grid on

nexttile; hold on
for timeValue = [1 2 3 5]
    [~,j] = min(abs(selected.time_s-timeValue));
    plot(wn,absorbance(j,:),'LineWidth',1.5, ...
        'DisplayName',sprintf('%.1f s',selected.time_s(j)));
end
xlabel('Wavenumber (cm^{-1})'); ylabel('Synthetic absorbance');
title('Acetone-only synthetic spectrum'); legend('Location','best'); grid on

exportgraphics(fig,fullfile(outputDir,'first_forward_curves_matlab.png'),'Resolution',220);

slowTable = table(slow.time_h(:),slow.r(:),slow.B_mM(:),slow.Q_mM(:), ...
    slow.A_blood_mM(:),'VariableNames', ...
    {'time_h','r_dimensionless','BHB_mM','AcAc_mM','blood_acetone_mM'});
writetable(slowTable,fullfile(outputDir,'slow_metabolism_matlab.csv'));

save(fullfile(outputDir,'forward_results.mat'),'p','slow','breaths', ...
    'breathHours','wn','absorbance');

results = struct('parameters',p,'slow',slow,'breaths',{breaths}, ...
    'wavenumber_cm1',wn,'absorbance',absorbance);
fprintf('Final BHB: %.3f mmol/L\n',slow.B_mM(end));
fprintf('Final blood acetone: %.1f micromol/L\n',1000*slow.A_blood_mM(end));
fprintf('Five-hour breath plateau: %.3f ppm\n',selected.C_exh_ppm(end));
end
