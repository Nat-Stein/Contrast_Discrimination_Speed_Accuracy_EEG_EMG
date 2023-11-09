%% Figure2_plot


% Set up layout of the figure
clear allWidth
topoWidth = [0 800]; SEMwidth = [0 300]; 
allWidth = {topoWidth, xlimSL, SEMwidth, xlimRL};
allHeight = {0.9 0.6 0.6 0.6};
fHeight = 0.15;
x1st = 0.02; xLst = 0.04; plotDist = 0.05;
plotWidth = 1-x1st - xLst - (size(allWidth,2)-1)*plotDist;
xWidths(1)=0;
for f=1:size(allWidth,2)
    xWidths(f+1) = allWidth{f}(2)-allWidth{f}(1);
end
totalX = sum(xWidths);
yPos = [0.08 0.8];
% ---------------------------------
t_bin = [250 450]; % Stats window for error bar plots (SEM)

%% Plot topography of SSVEP (average of 20Hz and 25Hz)

% Values are computed in Figure2_compute
load(fullfile(figData, 'SSVEPs'), 'SSVEPtopoSubtrBL')
load(fullfile(figData, 'stats_RL_SsvepPupil'))
load(fullfile(figData, 'allF_SSVEP_SL')) % Testing temporal extend of significant effects on SSVEP

figure('Position', [100, 100, 1000, 600]); hold on

yy=1;
% Define position and dimenstions of subplot
f=1; plotStart = x1st+plotWidth*sum(xWidths(1:f))/totalX + (f-1)*plotDist; plotX = plotWidth*xWidths(f+1)/totalX;
plotY = yPos(yy);
subplot('Position',[plotStart plotY plotX allHeight{f}]); hold on
% subplot(1,3,1); hold on
avgtopo = nan(17,97);
for subj = subjects
    trials =  goodTrialsComb{subj};
    avgtopo(subj,:) = squeeze(nanmean(nanmean(SSVEPtopoSubtrBL{subj}(1:2,:,trials),1),3));
end
topoplot(nanmean(avgtopo,1),chanlocs,'electrodes','off','shading','interp','style','straight','plotrad',0.56,'maplimits',[0 0.8]);
title(['SSVEP topography'])
axprefs(gca)

% ----------------------------------------------------------
% 
% Define position and dimenstions of subplot
f=2; plotStart = x1st+plotWidth*sum(xWidths(1:f))/totalX + (f-1)*plotDist; plotX = plotWidth*xWidths(f+1)/totalX;
plotY = yPos(yy);
subplot('Position',[plotStart plotY plotX allHeight{f}]); hold on

rw=1; 
% Stimulus-locked
% subplot(2,3,5); hold on
for cc = 1:2
    for l = 1:2
        for lr = 1:2
            plotSSVEP = squeeze(meanSL_SSVEP(cc,l,lr,:,rw)) -...
                squeeze(nanmean(nanmean(nanmean(meanSL_SSVEP(:,:,:,find(T==-fftlenSSVEP/2),rw),3),2),1));
            plot(T,plotSSVEP,'Color',colores{3-l,cc},'LineStyle',dash{1},'Linewidth',lr)
        end
    end
end
xlim(xlimSL);  line([0 0],[ylimitsSSVEP],'Color','k'); ylim(ylimitsSSVEP);
line(350*[1 1],ylimitsSSVEP,'Color','k');
line(550*[1 1],ylimitsSSVEP,'Color','k');
title(['Stimulus-locked SSVEP']); 
xlabel('Time [ms]'); ylabel('SSVEP power [uV]')
axprefs(gca)

% ---------------------------------------------
% Mean SSVEP amplitude in stats window with SEM across subjects

% Define position and dimenstions of subplot
f=3; plotStart = x1st+plotWidth*sum(xWidths(1:f))/totalX + (f-1)*plotDist; plotX = plotWidth*xWidths(f+1)/totalX;
plotY = yPos(yy);
subplot('Position',[plotStart plotY plotX allHeight{f}]); hold on

clear dataMatrix dataOrig subjMean
dataOrig = squeeze(nanmean(meanSL_SSVEPS(:,:,:,:,subjects,find(T>= t_bin(1) & T<= t_bin(2))),6)) -...
    squeeze(nanmean(nanmean(nanmean(meanSL_SSVEP(:,:,:,find(T==-fftlenSSVEP/2),rw),3),2),1)); % SSVEP at baseline - same for all conditions
subjMean = squeeze(nanmean(nanmean(nanmean(dataOrig,3),2),1));
for s = size(dataOrig,4)
    dataMatrix(:,:,:,s) = squeeze(dataOrig(:,:,:,s)) - squeeze(repmat(subjMean(s),2,2,2));
end
xpos = 0; 
for l=1:2
    xpos = xpos+0.2;
    for cc = 1:2
        xpos = xpos+0.2;
        for lr = 1:2
            v = nanstd(dataMatrix(cc,l,lr,:))./sqrt(size(dataOrig,4));
            m = nanmean(dataOrig(cc,l,lr,:),4);
            
            plot(xpos,m,'*','Color',colores{3-l,cc},'LineStyle',dash{1},'Linewidth',2)
            line(xpos*[1 1],[m-v m+v],'Color',colores{3-l,cc},'LineStyle',dash{1},'Linewidth',lr)
        end
    end
end
ylim(ylimitsSSVEP)
xlim([0 1.6])
title('250-450ms')
axprefs(gca)

% ----------------------------------------------------------
% Response-locked SSVEP traces
% Define position and dimenstions of subplot
f=4; plotStart = x1st+plotWidth*sum(xWidths(1:f))/totalX + (f-1)*plotDist; plotX = plotWidth*xWidths(f+1)/totalX;
plotY = yPos(yy);
subplot('Position',[plotStart plotY plotX allHeight{f}]); hold on
% Plot
for cc = 1:2
    for l = 1:2
        for lr = 1:2
            plotSSVEP = squeeze(nanmean(meanRL_SSVEPS(cc,l,lr,rw,:,:),5)) -...
                squeeze(nanmean(nanmean(nanmean(meanSL_SSVEP(:,:,:,find(T==-fftlenSSVEP/2),rw),3),2),1));
            plot(Tr,plotSSVEP,'Color',colores{3-l,cc},'LineStyle',dash{1},'Linewidth',lr)
        end
    end
end
xlim(xlimRL);  line([0 0],[ylimitsSSVEP],'Color','k'); ylim(ylimitsSSVEP);
title(['Response-locked SSVEP']);
xlabel('Time [ms]'); ylabel('SSVEP power [uV]')
axprefs(gca)

% ----------------------------------------------------------
% F-values that track significance of effect of speed pressure over time
% Define position and dimenstions of subplot
yy = 2; f=2; plotStart = x1st+plotWidth*sum(xWidths(1:f))/totalX + (f-1)*plotDist; plotX = plotWidth*xWidths(f+1)/totalX;
plotY = yPos(yy);
subplot('Position',[plotStart plotY plotX fHeight]); hold on
% Plot F-values of effects of Left/Right Target, Left/Right *
% Speed/Accuracy, and Left/Right * Contrast on SSVEP at every individual 
% time point 
cols = {[0, 0.5, 0], [0, 0.75, 0.3], [0.9290, 0.6940, 0.1250]};
plotEff = [5]; % plotEff = [3 5 6];
% Plot
for f = 1:length(plotEff)
plot(T,allF{plotEff(f)}(:,1),'Color',cols{f},'LineStyle',dash{1},'Linewidth',2)
end
line([T(1) T(end)],4.5431*[1 1],'Color','k');
line(350*[1 1],[0 15],'Color','k');
line(550*[1 1],[0 15],'Color','k');
legend('Left/Right * Speed/Accuracy')
% legend('Left/Right Target','Left/Right * Speed/Accuracy','Left/Right * Contrast')
xlim(xlimSL); ylabel('F-value'); 
axprefs(gca)

figure; hold on 
for f = 1:length(plotEff)
plot(T,allF{plotEff(f)}(:,2),'Color',cols{f},'LineStyle',dash{1},'Linewidth',2)
end
line([T(1) T(end)],0.05*[1 1],'Color','k');
line(350*[1 1],[0 0.1],'Color','k');
line(550*[1 1],[0 0.1],'Color','k');
legend('Left/Right Target','Left/Right * Speed/Accuracy','Left/Right * Contrast')
xlim(xlimSL); ylim([0 0.1]); ylabel('p-value'); xlabel('Time [ms]')
axprefs(gca)


%% Plot SSVEP for Correct vs Error trials
% Low Contrast only, because several subjects made very few or no mistakes
% on the high Contrast trials (especially for Accuracy condition)
l = 1;

load(fullfile(figData, 'SSVEP_RW'))

% Values are computed in Figure2_compute
% Stimulus-locked SSVEP for correct vs wrong trials
dash = {'-','--','-.'};
plt=0; clear legs
figure; hold on
for lr = 1:2
    for cc = 2:-1:1
        for rw = 1:2
            for p = 1:numbins
                plotSSVEP = squeeze(nanmean(meanSL_SSVEPS_bins(cc,l,lr,rw,:,:,p),5)) -...
                    squeeze(nanmean(nanmean(nanmean(nanmean(nanmean(nanmean(nanmean(meanSL_SSVEPS_bins(:,:,lr,:,:,find(T==-fftlenSSVEP/2),:),7),6),5),4),3),2),1));
                % Subtraction of SSVEP at baseline is only done for
                % plotting purposes
                plot(T,plotSSVEP,'Color',colores{3-l,cc},'LineStyle',dash{rw},'Linewidth',lr)
                plt=plt+1; legs{plt} = [CorrWrong{2,rw} ' ' LoHi{l} ' ' LeftRight{lr} ' ' AccSpd{cc}];
            end
        end
    end
end
legend(legs)
title('SSVEP correct vs wrong')
ylabel('SSVEP left-right')
xlabel('time')



