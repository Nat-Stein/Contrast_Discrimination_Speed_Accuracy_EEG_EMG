% SuppFigure6_plot

nbins = 3;

%% Plot traces for 3 RT bins and indicate significant correlations with stars

load(fullfile(figData, [num2str(nbins) 'RTbin_signs']))
load(fullfile(figData, ['CPP_MB+slope_' num2str(nbins) 'RTbin']))
% -----------------------------------------
% Define plotting parameters
c2=4;
hCPP = -3;
hMB = 9.0;
gCol = 0.4;
ylimSSVEP = [-1 2.5];
ylimCPP = [-5 25];
ylimMB = [6.5 9.5];
hSSVEP = -0.5;
smoothWidth = 3;

% -----------------------------------------
% Set up figure
clear allWidth xWidths
xlimRL = [-800 200]; xlimSL = [-200 1500]; widthBL = [0 500]; xlimRT = [400 1400];
allWidth = {xlimSL, xlimRL};
x1st = 0.05; xLst = 0.04; plotDist = 0.05;
plotWidth = 1-x1st - xLst - (size(allWidth,2)-1)*plotDist;
xWidths(1)=0;
for f=1:size(allWidth,2)
    xWidths(f+1) = allWidth{f}(2)-allWidth{f}(1);
end
totalX = sum(xWidths);

yPos = [0.08 0.40 0.72]; yHeight = 0.25;

figure('Position', [100, 100, 700, 800]); hold on

% -----------------------------------------
% Compute mean RTs
meanBinRTs = squeeze(nanmean(nanmean(nanmean(meanRTs,3),2),4));

% -----------------------------------------
% Stim-locked SSVEP
yy=1; f=1; plotStart = x1st+plotWidth*sum(xWidths(1:f))/totalX + (f-1)*plotDist; plotX = plotWidth*xWidths(f+1)/totalX;
plotY = yPos(yy);
subplot('Position',[plotStart plotY plotX yHeight]); hold on
for b = 1:nbins
    plot(T,smooth(squeeze(nanmean(nanmean(nanmean(bin2SSVEPslCs(b,:,:,:,:),2),3),4)),smoothWidth),...
        'Color',colores{4-b,c2},'LineStyle',dash{1},'Linewidth',2)
end
plot(T,(hSSVEP)*p_RT{yy,f}(1,:),'*','Color',gCol*[1 1 1]);

title(['SSVEP stimulus-locked, ' num2str(nbins) ' RT-bins'])
xlabel('Times [ms]'); ylabel('Amplitude')
xlimSL = [-200 1500];
xlim(xlimSL); ylim(ylimSSVEP)
line([0 0],ylimSSVEP,'Color','k')
for b = 1:nbins
    line(meanBinRTs(b)*1000*[1 1],ylimSSVEP,'Color',colores{4-b,c2})
end
legend('fast','medium','slow','Location','SouthEast');

% -----------------------------------------
% Response-locked SSVEP
f=2; plotStart = x1st+plotWidth*sum(xWidths(1:f))/totalX + (f-1)*plotDist; plotX = plotWidth*xWidths(f+1)/totalX;
plotY = yPos(yy);
subplot('Position',[plotStart plotY plotX yHeight]); hold on
for b = 1:nbins
    plot(Tr,smooth(squeeze(nanmean(nanmean(nanmean(bin2SSVEPrlCs(b,:,:,:,:),2),3),4)),smoothWidth),...
        'Color',colores{4-b,c2},'LineStyle',dash{1},'Linewidth',2)
end
plot(Tr,(hSSVEP)*p_RT{yy,f}(1,:),'*','Color',gCol*[1 1 1]);

title(['SSVEP response-locked, ' num2str(nbins) ' RT-bins'])
xlabel('Times [ms]'); ylabel('Amplitude')
xlimRL = [-800 200];
line([0 0],ylimSSVEP,'Color','k')
xlim(xlimRL); ylim(ylimSSVEP)

% -----------------------------------------
% Stim-locked CPP
yy=2; f=1; plotStart = x1st+plotWidth*sum(xWidths(1:f))/totalX + (f-1)*plotDist; plotX = plotWidth*xWidths(f+1)/totalX;
plotY = yPos(yy);
subplot('Position',[plotStart plotY plotX yHeight]); hold on
for b = 1:nbins
    plot(t,smooth(squeeze(nanmean(nanmean(nanmean(bin2CPPslCs(b,:,:,:,:),2),3),4)),smoothWidth),...
        'Color',colores{4-b,c2},'LineStyle',dash{1},'Linewidth',2)
end
plot(t,(hCPP)*p_RT{yy,f}(1,:),'*','Color',gCol*[1 1 1]);

title(['CPP stimulus-locked, ' num2str(nbins) ' RT-bins'])
xlabel('Times [ms]'); ylabel('Amplitude')
xlimSL = [-200 1500];
xlim(xlimSL); ylim(ylimCPP)
line([0 0],ylimitsCPP,'Color','k')
for b = 1:nbins
    line(meanBinRTs(b)*1000*[1 1],ylimCPP,'Color',colores{4-b,c2})
end
legend('fast','medium','slow');

% -----------------------------------------
% Response-locked CPP
f=2; plotStart = x1st+plotWidth*sum(xWidths(1:f))/totalX + (f-1)*plotDist; plotX = plotWidth*xWidths(f+1)/totalX;
plotY = yPos(yy); x = [-350:-50]; mid = (x(1)+x(end))/2;
subplot('Position',[plotStart plotY plotX yHeight]); hold on
for b = 1:nbins
    plot(tr,smooth(squeeze(nanmean(nanmean(nanmean(bin2CPPrlCs(b,:,:,:,:),2),3),4)),smoothWidth),...
        'Color',colores{4-b,c2},'LineStyle',dash{1},'Linewidth',2)
end
plot(tr,hCPP*p_RT{yy,f}(1,:),'*','Color',gCol*[1 1 1]);

title(['CPP response-locked,' num2str(nbins) ' RT-bins'])
xlabel('Times [ms]'); ylabel('Amplitude')
xlimRL = [-800 200];
line([0 0],ylimCPP,'Color','k')
xlim(xlimRL); ylim(ylimCPP)

% -----------------------------------------
% Stim-locked Mu/Beta
yy=3; f=1; plotStart = x1st+plotWidth*sum(xWidths(1:f))/totalX + (f-1)*plotDist; plotX = plotWidth*xWidths(f+1)/totalX;
plotY = yPos(yy);
subplot('Position',[plotStart plotY plotX yHeight]); hold on
set(gca,'Ydir','reverse');
for b = 1:nbins
    plot(T,smooth(squeeze(nanmean(nanmean(nanmean(bin2MBslCs(b,:,:,:,:),2),3),4)),smoothWidth),...
        'Color',colores{4-b,c2},'LineStyle',dash{1},'Linewidth',2)
    
    plot(T,smooth(squeeze(nanmean(nanmean(nanmean(bin2MBslCsI(b,:,:,:,:),2),3),4)),smoothWidth),...
        'Color',colores{4-b,c2},'LineStyle',dash{2},'Linewidth',2)
end
for ic = 1:2
    plot(T,(hMB+ic*0.1)*p_RT{yy,f}(2,:),'*','Color',gCol*[1 1 1]/ic);
end

title(['MB stimulus-locked,' num2str(nbins) ' RT-bins'])
xlabel('Times [ms]'); ylabel('Amplitude')
xlimSL = [-200 1500];
xlim(xlimSL); ylim(ylimMB)
line([0 0],ylimitsMB,'Color','k')
for b = 1:nbins
    line(meanBinRTs(b)*1000*[1 1],ylimMB,'Color',colores{4-b,c2})
end
legend('fast exec.','fast withh.','medium exec.','medium withh.',...
    'slow exec.','slow withh.',...
    'sign. exec.','sign. withh.',...
    'Location','NorthEast');

% -----------------------------------------
% Response-locked Mu/Beta
f=2; plotStart = x1st+plotWidth*sum(xWidths(1:f))/totalX + (f-1)*plotDist; plotX = plotWidth*xWidths(f+1)/totalX;
plotY = yPos(yy); x = [-350:50:-150]; mid = (x(1)+x(end))/2;
subplot('Position',[plotStart plotY plotX yHeight]); hold on
set(gca,'Ydir','reverse');
for b = 1:nbins
    plot(Tr,smooth(squeeze(nanmean(nanmean(nanmean(bin2MBrlCs(b,:,:,:,:),2),3),4)),smoothWidth),...
        'Color',colores{4-b,c2},'LineStyle',dash{1},'Linewidth',2)
    
    plot(Tr,smooth(squeeze(nanmean(nanmean(nanmean(bin2MBrlCsI(b,:,:,:,:),2),3),4)),smoothWidth),...
        'Color',colores{4-b,c2},'LineStyle',dash{2},'Linewidth',2)
    
end
title(['MB response-locked, ' num2str(nbins) ' RT-bins'])
xlabel('Times [ms]'); ylabel('Amplitude')
for ic = 1:2
    plot(Tr,(hMB+ic*0.1)*p_RT{yy,f}(ic,:),'*','Color',gCol*[1 1 1]/ic);
end
xlimRL = [-800 200];
line([0 0],ylimMB,'Color','k')
xlim(xlimRL); ylim(ylimMB)


