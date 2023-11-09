% Figure3_plot

load(fullfile(figData, 'stats_RL_SsvepPupil'))

plt = 0; clear legs
t_pup = t_targ(1:size(meanSL_PUPIL,4)); rw = 1;
figure; hold on
for cc = 1:2
    for l = 1:2
        plot(t_pup,squeeze(meanSL_PUPIL(cc,l,rw,:))./pupNorm,...
            'Color',colores{3-l,cc},'LineStyle',dash{1},'Linewidth',2)
        plt=plt+1; legs{plt} = [LoHi{l} ' ' AccSpd{cc}];
        
    end
end
xlim(xlimSL);
line([0 0],[pupilLim],'Color','k'); % ylim(pupilLim);
axprefs(gca)
title(['Stimulus-locked Pupil size pooled']);
xlabel('Time [ms]'); ylabel('Pupil area [a. u.]')
legend(legs,'Location','SouthEast')