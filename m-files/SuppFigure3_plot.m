% SuppFigure3_plot


%% Load matrices needed for plotting all traces
load(fullfile(figData, 'plot_SuppFig3_CPPaep'))

%%
% Set up layout of the figure
clear allWidth
allWidth = {widthBL, xlimSL, xlimRL,xlimRT};
x1st = 0.02; xLst = 0.04; plotDist = 0.05;
plotWidth = 1-x1st - xLst - (size(allWidth,2)-1)*plotDist;
xWidths(1)=0;
for f=1:size(allWidth,2)
    xWidths(f+1) = allWidth{f}(2)-allWidth{f}(1);
end
totalX = sum(xWidths);
yPos = 0.15 ; yHeight = 0.7;
% ---------------------------------


figure('Position', [100, 100, 1200, 400]); hold on

% --------------------------------------------------------------------
% Plot CPP topography (at response)
yy=1;
% Define position and dimenstions of subplot
f=1; plotStart = x1st+plotWidth*sum(xWidths(1:f))/totalX + (f-1)*plotDist; plotX = plotWidth*xWidths(f+1)/totalX;
plotY = yPos(yy);
subplot('Position',[plotStart plotY plotX yHeight]); hold on
% topoplot is a eeblab function
topoplot(nanmean(avgtopo,1),chanlocs,'maplimits',[-5 10],'shading','interp','style','straight','plotrad',0.56); title('CPP at RT')
axprefs(gca)

% --------------------------------------------------------------------
% Stimulus-locked CPP
% Define position and dimenstions of subplot
f=2; plotStart = x1st+plotWidth*sum(xWidths(1:f))/totalX + (f-1)*plotDist; plotX = plotWidth*xWidths(f+1)/totalX;
plotY = yPos(yy);
subplot('Position',[plotStart plotY plotX yHeight]); hold on
% Compute and plot traces for all four conditions
for cc = 1:2
    for l = 1:2
        plot(t,squeeze(meanSL_CPP(cc,l,:)),'Color',colores{3-l,cc},'LineStyle',dash{1},'Linewidth',2)
    end
end
axprefs(gca)
xlim(xlimSL); ylim(ylimitsCPP); line([0 0],[ylimitsCPP],'Color','k')
title('Stimulus-locked CPP'); xlabel('Time [ms]'); ylabel('CPP [uV]')

% --------------------------------------------------------------------
% Response-locked CPP
% Define position and dimenstions of subplot
f=3; plotStart = x1st+plotWidth*sum(xWidths(1:f))/totalX + (f-1)*plotDist; plotX = plotWidth*xWidths(f+1)/totalX;
plotY = yPos(yy);
subplot('Position',[plotStart plotY plotX yHeight]); hold on
% Compute and plot traces for all four conditions
for cc = 1:2
    for l = 1:2
        plot(tr,squeeze(meanRL_CPP(cc,l,:)),'Color',colores{3-l,cc},'LineStyle',dash{1},'Linewidth',2)
    end
end
line([t_RT(1) t_RT(2)],[ylimitsCPP(2)-2 ylimitsCPP(2)-2],'Color',0.3*[1 1 1],'Linewidth',2)
axprefs(gca); xlim(xlimRL); ylim(ylimitsCPP);
line([0 0],[ylimitsCPP],'Color','k')
title('Response-locked CPP'); xlabel('Time [ms]'); ylabel('CPP [uV]')


% --------------------------------------------------------------------
% CPP at decision commitment (-130 to -70ms before click) plotted over RT
% Define position and dimenstions of subplot
f=4; plotStart = x1st+plotWidth*sum(xWidths(1:f))/totalX + (f-1)*plotDist; plotX = plotWidth*xWidths(f+1)/totalX;
plotY = yPos(yy);
subplot('Position',[plotStart plotY plotX yHeight]); hold on


% Plot standard error of the mean (across subjects) for all conditions
% and RT bins
clear dataMatrix dataMatrixRT
subjMean = squeeze(nanmean(nanmean(nanmean(nanmean(nanmean(meanCPPCtr,6),5),3),2),1));
subjMeanRT = squeeze(nanmean(nanmean(nanmean(nanmean(nanmean(meanRTCtr,6),5),3),2),1));
for subj = subjects
    dataMatrix(:,:,:,subj) = squeeze(nanmean(nanmean(meanCPPCtr(:,:,:,subj,:,:),6),5)) - squeeze(repmat(subjMean(subj),CPPnumbins,2,2));
    dataMatrixRT(:,:,:,subj) = squeeze(nanmean(nanmean(meanRTCtr(:,:,:,subj,:,:),6),5)) - squeeze(repmat(subjMeanRT(subj),CPPnumbins,2,2));
end

for cc=1:2
    for l = 1:2
        % Plot averages within conditions as a line
        plot(squeeze(nanmean(nanmean(nanmean(meanRTCtr(:,cc,l,:,:,:),6),5),4)),...
            squeeze(nanmean(nanmean(nanmean(meanCPPCtr(:,cc,l,:,:,:),6),5),4)),...
            'Color',colores{3-l,cc},'LineStyle',dash{1},'Linewidth',2)
        
        % Plot SEM in RT and signal amplitude as lines
        for b=1:CPPnumbins
            v = nanstd(dataMatrix(b,cc,l,:))./sqrt(length(subjects));
            m = nanmean(nanmean(nanmean(nanmean(meanCPPCtr(b,cc,l,:,:,:),6),5),4));
            
            vRT = nanstd(dataMatrixRT(b,cc,l,:))./sqrt(length(subjects));
            mRT = nanmean(nanmean(nanmean(nanmean(meanRTCtr(b,cc,l,:,:,:),6),5),4));
            
            line([mRT mRT],[m-v m+v],'Color',colores{3-l,cc},'LineStyle',dash{1},'Linewidth',1)
            
            
            line([mRT-vRT mRT+vRT],[m m],'Color',colores{3-l,cc},'LineStyle',dash{1},'Linewidth',1)
        end
    end
end
axprefs(gca); xlim(xlimRT); ylim(ylimitsCPP)
title(['CPP over RT (' num2str(t_RT(1)) ' to ' num2str(t_RT(2)) 'ms)']); xlabel('Time [ms]'); ylabel('CPP [uV]')

