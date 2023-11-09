% SuppFigure2_plot

%% Plot single-frequency SSVEP and the F-values of the effect of 
% Speed/Accuracy condition for 20, 22.5 and 25Hz

load(fullfile(figData, 'SSVEPindFreq'))
load(fullfile(figData, 'SSVEPindFreq_F'),'allPs','allFs')


plotF = [20 22.5 25]; ylimitsFreqs = [4 13];
figure; hold on; plt=length(plotF);
for freq = 1:length(plotF)
    plt=plt+1;subplot(2,length(plotF),plt); hold on
    for l = 1:2
        for cc = 1:2
            for lr = 1:2
                
                plotSSVEP = squeeze(nanmean(meanSL_FreqSubj(cc,l,lr,find(Fssvep==plotF(freq)),:,:),6));
                plot(T,plotSSVEP,'Color',colores{3-l,cc},'LineStyle',dash{1},'Linewidth',lr)
                
            end
            mRT = nanmean(nanmean(RTmeanS(cc,l,:,:),3),4);
            line(mRT*1000*[1 1],ylimitsFreqs,'Color',colores{3-l,cc},'LineStyle',dash{2},'Linewidth',1.5)
        end
    end
    axis tight; ylim(ylimitsFreqs); line([0 0],[ylimitsFreqs],'Color','k'); xlabel('Time [ms]'); ylabel('SSVEP power [uV]')
    axprefs(gca);
end
% legend('Left, Acc, lowC','Left, Acc, highC','Left, Spd, lowC','Left, Spd, highC',...
%     'Right, Acc, lowC','Right, Acc, highC','Right, Spd, lowC','Right, Spd, highC')
% Plot F-values
take_varis = [1]; thisLim = [0 15];
% take_varis = [1 3 6]; thisLim = [0 30];
for freq = 1:length(take_freqs)
   subplot(2,3,freq); hold on 
   plot(T,squeeze(allFs(freq,:,take_varis)),'LineWidth',2)
   line([T(1) T(end)],[4.45 4.45],'Color','k')
   line([0 0],thisLim,'Color','k')
   axis tight; title([num2str(take_freqs(freq)) 'Hz'])
   xlabel('Time [ms]'); ylabel('F-value')
   axprefs(gca); 
end
legend(output.labels{take_varis})

