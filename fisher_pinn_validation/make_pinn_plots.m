function make_pinn_plots(P,outDir,mode)
h=P.temporal_holdout;
fig=figure('Color','w','Position',[80 80 1050 430]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; scatter(P.BHB_true_mM(h),P.BHB_pred_mM(h),35,'filled'); hold on
lims=[min([P.BHB_true_mM(h);P.BHB_pred_mM(h)]) max([P.BHB_true_mM(h);P.BHB_pred_mM(h)])];
plot(lims,lims,'k--'); axis equal; xlim(lims); ylim(lims); grid on
xlabel('True BHB (mmol/L)'); ylabel('PINN BHB (mmol/L)'); title(mode+" Gamma: holdout");
nexttile; scatter(P.slope_true_mM_h(h),P.slope_pred_mM_h(h),35,'filled'); hold on
lims=[min([P.slope_true_mM_h(h);P.slope_pred_mM_h(h)]) max([P.slope_true_mM_h(h);P.slope_pred_mM_h(h)])];
plot(lims,lims,'k--'); axis equal; xlim(lims); ylim(lims); grid on
xlabel('True dB/dt'); ylabel('PINN dB/dt'); title('Trajectory derivative');
exportgraphics(fig,fullfile(outDir,'pinn_holdout.png'),'Resolution',220);

ids=unique(P.patient_id); fig2=figure('Color','w','Position',[80 80 1100 650]);
tiledlayout(2,ceil(numel(ids)/2),'TileSpacing','compact','Padding','compact');
for k=1:numel(ids)
 Q=P(P.patient_id==ids(k),:); nexttile;
 plot(Q.time_h,Q.BHB_true_mM,'k-','LineWidth',1.5); hold on
 plot(Q.time_h,Q.BHB_pred_mM,'b-','LineWidth',1.5);
 scatter(Q.time_h(Q.temporal_holdout),Q.BHB_pred_mM(Q.temporal_holdout),38,'filled');
 grid on; title(sprintf('Patient %d: %s',ids(k),string(Q.scenario(1))));
 xlabel('Time (h)'); ylabel('BHB (mmol/L)');
end
legend('Truth','PINN','Holdout','Location','best');
exportgraphics(fig2,fullfile(outDir,'pinn_trajectories.png'),'Resolution',220);
end
