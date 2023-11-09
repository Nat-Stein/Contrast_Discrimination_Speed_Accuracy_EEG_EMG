% Figure6_plot

% ----------------------------------------------
% Plot layout
clear allWidth xWidths
xlimRL = [-800 200]; xlimSL = [-200 1500]; widthBL = [0 1000]; xlimRT = [400 1400];
allWidth = {widthBL, xlimRL, xlimRT};
x1st = 0.05; xLst = 0.04; plotDist = 0.05;
plotWidth = 1-x1st - xLst - (size(allWidth,2)-1)*plotDist;
xWidths(1)=0;
for f=1:size(allWidth,2)
    xWidths(f+1) = allWidth{f}(2)-allWidth{f}(1);
end
totalX = sum(xWidths);
yPos = [0.08 0.58]; yHeight = 0.38;


%% EMG activation for 2 RT-bins over RT (different RT bins)
% Distribution of EMB onset bursts, response-locked EMG amplitude of 2 RT-bins
% and EMG amplitude at response plotted for 4 RT bins

load(fullfile(figData, 'EMG_plotMat'))

% When to have the middle of the 'at response' window (width = 100ms). Steps of 25ms!!!
figure; hold on; plt=0; clear legs
numbins = 2; percentiles = [linspace(0,100,numbins+1)]; numbinsRT= 4; 
for cho = 1:2
    
    %------------------------------
    % Distribution of EMG onset bursts over the course of a trial
    % Define plot layout
    yy=3-cho;
    f=1; plotStart = x1st+plotWidth*sum(xWidths(1:f))/totalX + (f-1)*plotDist; plotX = plotWidth*xWidths(f+1)/totalX;
    plotY = yPos(yy);
    subplot('Position',[plotStart plotY plotX yHeight]); hold on
    
    for cc = 1:2
        thisBin = tr([find(tr==histOnsets(1)):binWidthEMG(cho):find(tr==histOnsets(2))]);
        
        % Plot traces
        plot(-thisBin,squeeze(emgOnsetStats(cho,cc,:)),'Color',colores{1,cc},'LineStyle',dash{1},'Linewidth',2);
        plt = plt+1; legs{plt} = [AccSpd{cc}];
        title(['EMG onsets ' choName{cho}]); xlabel(['Time [ms], binWidthEMG = ' num2str(binWidthEMG(cho))]);
        ylabel('Mean number EMG onsets'); xlim([min(thisBin) max(thisBin)])
        legend(legs); xlim([-400 0]); ylim(ylimsOnsets{cho})
        
    end
    
    %------------------------------
    % EMG amplitude over the course of a trial
    % Define plot layout
    f=2; plotStart = x1st+plotWidth*sum(xWidths(1:f))/totalX + (f-1)*plotDist; plotX = plotWidth*xWidths(f+1)/totalX;
    plotY = yPos(yy);
    subplot('Position',[plotStart plotY plotX yHeight]); hold on
    
    for cc = 1:2
        
        % Plot traces
        for p = 1:length(percentiles)-1
            
            plot(Tr2,squeeze(EMG_SL(cho,cc,p,:)),'Color',colores{1,cc},'LineStyle',dash{1},'Linewidth',3-p); axis tight;
            
        end
    end
    % Plot description
    title(['Resp-locked ' choName{cho}]);
    xlabel('Time [ms]'); ylabel('Mean Power 10-250Hz')
    ylim(ylimitsEMG{cho}); xlim(xlimRL)
    line([0 0],ylimitsEMG{cho},'Color','k')
    
    %------------------------------
    % EMG amplitude at response plotted over RT
    % Define plot layout
    f=3; plotStart = x1st+plotWidth*sum(xWidths(1:f))/totalX + (f-1)*plotDist; plotX = plotWidth*xWidths(f+1)/totalX;
    plotY = yPos(yy);
    subplot('Position',[plotStart plotY plotX yHeight]); hold on
    
    % Compute SEM
    clear dataMatrix dataMatrixRT
    subjMean = squeeze(nanmean(nanmean(nanmean(nanmean(nanmean(meanEMG_4(:,cho,:,:,:),6),5),4),2),1));
    subjMeanRT = squeeze(nanmean(nanmean(nanmean(nanmean(meanRT_4,5),4),3),1));
    for subj = subjects;
        dataMatrix(:,:,subj) = squeeze(nanmean(meanEMG_4(:,cho,subj,:,:),1)) - squeeze(repmat(subjMean(subj),numbinsRT,2));
        dataMatrixRT(:,:,subj) = squeeze(meanRT_4(cho,subj,:,:)) - squeeze(repmat(subjMeanRT(subj),numbinsRT,2));
    end
    % Plot
    for cc=1:2
        plot(squeeze(nanmean(meanRT_4(cho,:,:,cc),2)*1000),squeeze(nanmean(nanmean(meanEMG_4(:,cho,:,:,cc),3),1)),...
            'Color',colores{3-l,cc},'LineStyle',dash{1},'Linewidth',2)
        
        for b=1:numbinsRT
            v = nanstd(dataMatrix(b,cc,:))./sqrt(length(subjects));
            m = nanmean(nanmean(meanEMG_4(:,cho,:,b,cc),3),1);
            
            vRT = nanstd(dataMatrixRT(b,cc,:)*1000)./sqrt(length(subjects));
            mRT = nanmean(nanmean(meanRT_4(cho,:,b,cc),2))*1000;
            
            line([mRT mRT],[m-v m+v],'Color',colores{3-l,cc},'LineStyle',dash{1},'Linewidth',1)
            line([mRT-vRT mRT+vRT],[m m],'Color',colores{3-l,cc},'LineStyle',dash{1},'Linewidth',1)
        end
        
    end
    
    title(['Amplitude at RT (' num2str(t_EMG{cho}) ') ' choName{cho}]); xlabel('Reaction time [ms]'); ylabel('Mean Power 10-250Hz')
    ylim(ylimitsEMG{cho}); xlim([400 1400])
    
end








