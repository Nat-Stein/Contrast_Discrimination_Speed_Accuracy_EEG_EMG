% SuppFigure1_plot


%% Points awarded for correct and incorrect choices as a function of RT 
% for all Speed/Accuracy conditions

figure('Position', [100, 100, 900, 900]); hold on
subplot(3,2,1); hold on
line([0 2.4],[60 60],'Color','b', 'LineWidth', 2)
line([0 2.4],[-60 -60],'Color','b','LineStyle','--', 'LineWidth', 2)
line([0 1],[100 100],'Color','r', 'LineWidth', 2)
line([1 2.4],[-120 -120],'Color','r', 'LineWidth', 2)
line([0 1],[-20 -20],'Color','r','LineStyle','--', 'LineWidth', 2)
xlim([0 2.4]); ylim([-140 140])
ylabel('Points'); xlabel('Reaction Time [s]'); title('Points per trial: Deadline conditions')

subplot(3,2,2); hold on
line([0 2.4],[70 60],'Color','b', 'LineWidth', 2)
line([0 2.4],[-70 -80],'Color','b','LineStyle','--', 'LineWidth', 2)
line([0 2.4],[100 -20],'Color','r', 'LineWidth', 2)
line([0 2.4],[-20 -120],'Color','r','LineStyle','--', 'LineWidth', 2)
xlim([0 2.4]); ylim([-140 140])
ylabel('Points'); xlabel('Reaction Time [s]'); title('Points per trial: Slope conditions')
legend('Accuracy correct', 'Accuracy error', 'Speed correct', 'Speed error')
%% Reaction time distribution for Deadline and Slope conditions separately
allHists = NaN*ones(4,2,2,length(histbins));

clear legs; plt = 0; % Prepare legend
for c = [1 3 2 4] % Four Speed/Accuracy conditions:
    % 1 = Accuracy Deadline, 2 = Speed Slope, 3 = Speed Deadline, 4 = Speed Slope
    
    for l = 1:2 % 1 = low contrast, 2 = high contrast
        this = NaN*ones(3,2,17,2,length(histbins));
        for d = 1:3 % Three delays (800, 1200 and 1600ms) between reward cue and evidence onset
            for lr = 1:2 % 1 = left, 2 = right
                % Devided by onset delays and left/right trials to then average
                % across both
                for subj = subjects
                    
                    % Correct trials of Speed/Accuracy condition c,
                    % contrast l, onset delay d, left/right lr
                    condsBeh{c,l,d,subj,1,lr} = find(indicators.onsedelay{subj} == d &...
                        indicators.cond{subj} == c &...
                        indicators.respLR{subj} == indicators.LR{subj} &...
                        indicators.ContrLevels{subj} == l &...
                        indicators.LR{subj} == lr  &...
                        indicators.RT{subj} >= minRT);
                    
                    % Incorrect trials of Speed/Accuracy condition c,
                    % contrast l, onset delay d, left/right lr
                    condsBeh{c,l,d,subj,2,lr} = find(indicators.onsedelay{subj} == d &...
                        indicators.cond{subj} == c &...
                        indicators.respLR{subj} ~= indicators.LR{subj} &...
                        indicators.ContrLevels{subj} == l &...
                        indicators.LR{subj} == lr  &...
                        indicators.RT{subj} >= minRT);
                    
                    for cw = 1:2 % 1 = repsonse correct, 2 = response incorrect
                        
                        % Proportion of trials in RT bin of each condition
                        this(d,lr,subj,cw,:) = hist(indicators.RT{subj}(condsBeh{c,l,d,subj,cw,lr}),histbins)./...
                            length([condsBeh{c,l,d,subj,2,lr} condsBeh{c,l,d,subj,1,lr}]);
                    end
                end
            end
        end
        for cw = 1:2 % 1 = repsonse correct, 2 = response incorrect
            
            subplot(3,2,DlSlo(c)+2); hold on; % Subplot 1 = Dedline conditions, Subplot 2 = Slope conditions
            
            allHists(c,l,cw,:) = squeeze(nanmean(nanmean(nanmean(this(:,:,:,cw,:),3),2),1));
            plot(histbins,squeeze(allHists(c,l,cw,:)),...
                'Color',colores{3-l,satCond(c)},'LineStyle',dash{cw},'LineWidth',3)
            
            plt = plt + 1; legs{plt} = [AccSpd{satCond(c)} ', ' LoHi{l} ', ' CorrWrong{2,cw}]; % Legend
        end
    end
    ylim([0 0.3]); xlim([0 2.4])
    ylabel('Proportion of trials');
    xlabel('Reaction Time [s]');
    title(['Reaction time distribution: ' DlSlo_name(DlSlo(c)) '  conditions'])
end
legend(legs)

%%

cc = 0; percentiles = [linspace(0,100,4+1)];
clear legs; plt = 0; % Prepare legend
for c = [1 3 2 4]
    for l = 1:2
        meanRT = NaN*ones(3,2,17,length(percentiles)); 
        meanCorr = NaN*ones(3,2,17,length(percentiles));
        for d = 1:3
            % Devided by onset delay to then average across onset delay bins
            for lr = 1:2
                for subj = subjects
                    
                    thisTrials = [condsBeh{c,l,d,subj,1,lr} condsBeh{c,l,d,subj,2,lr}]; numCorr = length(condsBeh{c,l,d,subj,1,lr});
                    thisRT = indicators.RT{subj}(thisTrials);
                    pc = [prctile(thisRT,percentiles)]; clear ts
                    for p = 1:length(pc)-1; ts{p} = find(thisRT>=pc(p) & thisRT<=pc(p+1)); meanRT(d,lr,subj,p) = nanmean(thisRT(ts{p}));
                        meanCorr(d,lr,subj,p) = length(find(ts{p}<=numCorr))./length(ts{p}); end
                end
            end
        end
        subplot(3,2,DlSlo(c)+4); hold on
        plot(squeeze(nanmean(nanmean(nanmean(meanRT,3),2),1)),squeeze(nanmean(nanmean(nanmean(meanCorr,3),2),1)),...
            'Color',colores{3-l,satCond(c)},'LineStyle',dash{1},'LineWidth',3)
        plt = plt + 1; legs{plt} = [AccSpd{satCond(c)} ', ' LoHi{l}]; % Legend
        
    end
    
    xlim([0 2.4]); ylim([0.5 1])
    ylabel('Proportion of correct trials'); xlabel('Reaction Time [s]'); title(['Conditional Accuracy: ' DlSlo_name{DlSlo(c)} ' conditions'])
end
legend(legs)

%% Compute Komolgorov-Smirnov to test whether RT distributions within 
% behavioural conditions are significantly different between Deadline
% and Slope conditions
for cc = 1:2
    for l = 1:2
        for cw = 1:2
            [h,p,ks2stat] = kstest2(squeeze(allHists(sats{cc}(1),l,cw,:)),squeeze(allHists(sats{cc}(2),l,cw,:)));
            disp([AccSpd{cc} ', ' LoHi{l} ', ' CorrWrong{2,cw} ' - h(Different distribution)=' num2str(h) ', p=' num2str(p)])
        end
    end
end

% Accuracy, Low Contrast, correct - h(Different distribution)=0, p=0.99503
% Accuracy, Low Contrast, wrong - h(Different distribution)=0, p=0.99503
% Accuracy, High Contrast, correct - h(Different distribution)=0, p=0.99503
% Accuracy, High Contrast, wrong - h(Different distribution)=0, p=0.71982
% Speed, Low Contrast, correct - h(Different distribution)=0, p=1
% Speed, Low Contrast, wrong - h(Different distribution)=0, p=0.71982
% Speed, High Contrast, correct - h(Different distribution)=0, p=0.91681
% Speed, High Contrast, wrong - h(Different distribution)=0, p=1







