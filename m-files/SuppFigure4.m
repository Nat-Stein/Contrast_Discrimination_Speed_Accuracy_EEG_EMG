% SuppFigure4

% ------------------------------------------------------------------
%% Test robustness of effect of RT and Speed/Accuracy emphasis on CPP 
% amplitude at response for different time windows of measurement

load(fullfile(figData, 'CPP_at_RT_twin'),'p_FE','chi2_FE','t_stats','wBin')
%% Compute response-locked unilateral motor potentials


load(fullfile(figData, 'meanRL_LRP'),'meanRL_LRP')

%% Plot figure

figure; hold on
subplot(2,1,2); hold on; rw=1;
for ic = 1:2
    plot(tr,squeeze(nanmean(nanmean(meanRL_LRP(:,:,ic,:,rw),2),1)),'Color',0.5*[1 1 1],'LineStyle',dash{ic},'Linewidth',2)
end
plot(tr,squeeze(nanmean(nanmean(meanRL_CPP,2),1)),'Color',0.1*[1 1 1],'LineStyle',dash{1},'Linewidth',2)
yl = [min(ylimitsLRP(1), ylimitsCPP(1)) max(ylimitsLRP(2), ylimitsCPP(2))];
xlim([-600 200]); ylim(yl); 
line([0 0],yl,'Color','k')
xlabel('Time [ms]'); ylabel('CPP, LRP [uV]')
legend('contralateral LRP','ipsilateral LRP','CPP','Location','NorthWest')
axprefs(gca)

% load(fullfile(dataFolder, conditions{1}, 'CPP_at_RT_twin'))
subplot(4,1,2); hold on;
plot(t_stats,p_FE,'LineWidth',2)
line([t_stats(1) t_stats(end)], 0.05*[1 1], 'Color', 0.3 * [1 1 1])
legend('RT','Speed/Accuracy','Location','NorthWest')
xlim([t_stats(1) t_stats(end)])
ylim([0 0.2])
axprefs(gca)
ylabel('p-value')


