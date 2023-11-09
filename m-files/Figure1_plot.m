% Figure1_plot



%% ------------------------------------------------------------
% Figure 1b - Response time distributions combined across Deadline and Slope
% conditions within behavioural condition
figure; hold on
subplot(2,1,1); hold on; clear legs; plt=0;
for cc= 1:2
    for l = 1:2
        this = NaN*ones(3,2,17,2,length(histbins));
        for d = 1:3
            % Devided by onset delay to then average across onset delay bins
            for lr = 1:2
                for subj = subjects
                    % name changed from condsBeh
                    condsBehPooled{cc,l,d,subj,1,lr} = find(indicators.onsedelay{subj} == d &...
                        (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                        indicators.respLR{subj}==indicators.LR{subj} &...
                        indicators.ContrLevels{subj} == l &...
                        indicators.LR{subj} == lr  &...
                        indicators.RT{subj}>=minRT);
                    
                    condsBehPooled{cc,l,d,subj,2,lr} = find(indicators.onsedelay{subj} == d &...
                        (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                        indicators.respLR{subj}~=indicators.LR{subj} &...
                        indicators.ContrLevels{subj} == l &...
                        indicators.LR{subj} == lr  &...
                        indicators.RT{subj}>=minRT);
                    
                    for cw = 1:2
                        % Plot RT histogram for every condition for correct
                        % and wrong trials
                        this(d,lr,subj,cw,:) = hist(indicators.RT{subj}(condsBehPooled{cc,l,d,subj,cw,lr}),histbins)./...
                            length([condsBehPooled{cc,l,d,subj,2,lr} condsBehPooled{cc,l,d,subj,1,lr}]);
                        
                    end
                end
            end
        end
        for cw = 1:2
            plot(histbins,squeeze(nanmean(nanmean(nanmean(this(:,:,:,cw,:),3),2),1)),...
                'Color',colores{3-l,cc},'LineStyle',dash{cw},'LineWidth',3)
            plt = plt + 1; legs{plt} = [AccSpd{cc} ', ' LoHi{l} ', ' CorrWrong{2,cw}]; % Legend
        end
    end
end
legend(legs)
ylim([0 0.3]); xlim([0 2.4])
ylabel('Proportion of trials'); xlabel('Reaction Time [s]'); title('Reaction time distribution: Deadline & slope pooled')



plotRTerror = 1; plotLine = 1;
clear condsBehCA
useBins = 10; useWins = [15]; useBinsP = [10];

subplot(2,1,2); hold on
percentiles = [linspace(0,100,useBins+1)];
meanRT = NaN*ones(2,17,length(percentiles),2,2); 
meanCorr = NaN*ones(2,17,length(percentiles),2,2);
for cc=1:2
    for l = 1:2
        % Pool across onset delays for lower within-condition noise
        for lr = 1:2
            for subj = subjects
                
                condsBehPooledCA{cc,l,subj,1,lr} = find((indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                    indicators.respLR{subj}==indicators.LR{subj} &...
                    indicators.ContrLevels{subj} == l &...
                    indicators.LR{subj} == lr &...
                    indicators.RT{subj}>=minRT);
                condsBehPooledCA{cc,l,subj,2,lr} = find((indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                    indicators.respLR{subj}~=indicators.LR{subj} &...
                    indicators.ContrLevels{subj} == l &...
                    indicators.LR{subj} == lr &...
                    indicators.RT{subj}>=minRT);
                
                % Separate the trials of each condition into RT bins
                thisTrials = [condsBehPooledCA{cc,l,subj,1,lr} condsBehPooledCA{cc,l,subj,2,lr}]; numCorr = length(condsBehPooledCA{cc,l,subj,1,lr});
                thisRT = indicators.RT{subj}(thisTrials);
                pc = [prctile(thisRT,percentiles)]; clear ts
                for p = 1:length(pc)-1; ts{p} = find(thisRT>=pc(p) & thisRT<=pc(p+1));
                    meanRT(lr,subj,p,cc,l) = nanmean(thisRT(ts{p})); 
                    meanCorr(lr,subj,p,cc,l) = length(find(ts{p}<=numCorr))./length(ts{p});
                end
            end
        end
        if plotLine == 1
            plot(squeeze(nanmean(nanmean(meanRT(:,:,:,cc,l),2),1)),squeeze(nanmean(nanmean(meanCorr(:,:,:,cc,l),2),1)),...
                'Color',colores{3-l,cc},'LineStyle',dash{1},'LineWidth',3)
        end
    end
end

% Add SEM error bars for each condition and RT bin
clear dataMatrix dataMatrixOns
subjMean = squeeze(nanmean(nanmean(nanmean(nanmean(meanCorr,5),4),3),1));
subjMeanRT = squeeze(nanmean(nanmean(nanmean(nanmean(meanRT,5),4),3),1));
dataMatrix = NaN*ones(size(meanCorr)); dataMatrixRT = NaN*ones(size(meanRT));
for subj = subjects;
    dataMatrix(:,subj,:,:,:) = squeeze(meanCorr(:,subj,:,:,:)) - squeeze(repmat(subjMean(subj),2,length(percentiles),2,2));
    dataMatrixRT(:,subj,:,:,:) = squeeze(meanRT(:,subj,:,:,:)) - squeeze(repmat(subjMeanRT(subj),2,length(percentiles),2,2));
end
for cc=1:length(sats)
    for l = 1:2
        for p = 1:length(pc)-1;
            v = nanstd(nanmean(dataMatrix(:,:,p,cc,l),1),[],2)./sqrt(length(subjects)); m = nanmean(nanmean(meanCorr(:,:,p,cc,l),1),2);
            mRT = nanmean(nanmean(meanRT(:,:,p,cc,l),1),2);
            line([mRT mRT],[m-v m+v],'Color',colores{3-l,cc},'LineStyle',dash{1},'Linewidth',2)
            
            if plotRTerror == 1
                v = nanstd(nanmean(dataMatrixRT(:,:,p,cc,l),1),[],2)./sqrt(length(subjects)); mRT = nanmean(nanmean(meanRT(:,:,p,cc,l),1),2);
                m = nanmean(nanmean(meanCorr(:,:,p,cc,l),1),2);
                line([mRT-v mRT+v],[m m],'Color',colores{3-l,cc},'LineStyle',dash{1},'Linewidth',2)
            end
        end
    end
end
xlim([0 2.4]); ylim([0.5 1])
ylabel('Proportion of correct trials'); xlabel('Reaction Time [s]'); 
title(['Conditional Accuracy: DL & SLP pooled - ' num2str(useBins) ' bins'])





