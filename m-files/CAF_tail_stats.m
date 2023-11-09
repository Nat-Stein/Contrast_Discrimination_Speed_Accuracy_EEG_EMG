%% CAF_tail_stats
% Measure slope of descending conditional accuracy past the point of peak
% performance to determine whether the drop in performance is greater under
% speed pressure

% To estimate each subject's time of peak accuracy, we compued the mean 
% accuracy in 6 RT percentiles (trials pooled across experimental 
% conditions). The time of peak accuracy was defined as the mean RT in the 
% RT-bin with the greatest response accuracy. 
% For each subject, we then fit lines to the "late" conditional accuracy
% functions within Speed/Accuracy emphasis, Left/Right trials and Contrast 
% level and tested whether Speed/Accuracy had a significant effect on the
% decrease in response accuracy over RT.

% To ensure that the results are robust to the number of RT-bins the time 
% of peak performance was computed based upon, we repeated the analysis 
% for 4 - 10 RT-bins

% Number of RT-bin the time of peak performance will be calculated for
useBins = [4 5 6 7 8 9 10]; 

allBestRTcut = nan(max(subjects),length(useBins));
minRT=0;

for binning = 1 : length(useBins)
    
    % Trials will be split into 4-10 [=useBins(binning)] percentiles based
    % on their RT
    percentiles = [linspace(0,100,useBins(binning)+1)];
    bestRTcut = nan(max(subjects),1);
    meanRTpool = nan(max(subjects),length(useBins));
    meanCorrPool = nan(max(subjects),length(useBins));
    
    % Determine best RT cut-off time per subject (and # of RT bins)
    for subj = subjects
        
        % Correct- and incorrect-response trials pooled within subject
        right = find(indicators.respLR{subj}==indicators.LR{subj} &...
            indicators.RT{subj}>=minRT);
        wrong = find(indicators.respLR{subj}~=indicators.LR{subj} &...
            indicators.RT{subj}>=minRT);
        
        % Sort trials into RT bins
        theseTrials = [right wrong]; 
        numCorr = length(right);
        thisRT = indicators.RT{subj}(theseTrials);
        pc = [prctile(thisRT,percentiles)]; clear ts
        for p = 1:length(pc)-1
            ts{p} = find(thisRT >= pc(p) & thisRT <= pc(p+1));
            meanRTpool(subj,p) = nanmean(thisRT(ts{p}));
            meanCorrPool(subj,p) = length(find(ts{p}<=numCorr))./length(ts{p});
        end
        
        mpRT = find(meanCorrPool(subj,:) == max(meanCorrPool(subj,:)));
        bestRT(subj) = mpRT(end);
        
        % bestRTcut specifies for each subject the mean RT of the RT-bin
        % with the peak performance
        bestRTcut(subj,:) = meanRTpool(subj,bestRT(subj));
        
    end
    
    allBestRTcut(:,binning) = bestRTcut;
    
    
    % Compute conditional accuracy for all trials with RTs slower than
    % cut-off RT
    cutRTs = repmat(allBestRTcut(:,binning),1,2,2); % Subject x Speed/Accuracy x Contrast
    disp('-------------------------------------------------')
    disp('3-Way repeated-measures ANOVA on descencing slope of CAF past peak performance')
    disp(['Cut-off RT based on ' num2str(useBins(binning)) ' RT bins'])
    
    % CAF_slope_stats splits trials with RT > cutRTs into 8 RT bins, fits
    % slopes to the data points within experimental conditions, computes
    % 3-Way ANOVA on the slopes and prints the results
    [lateCAF_anova] = CAF_slope_stats(indicators, cutRTs, subjects);
    
    if binning == 3
        % Save data for Bayes Factor analysis in JASP
        regression_table_CAF = [[1:size(lateCAF_anova,2)]; lateCAF_anova];
        csvwrite(fullfile(dataFolder, 'regression_table_CAF.csv'),regression_table_CAF)
        % Format:
        % Acc, lowC, Left - Acc, lowC, right - Acc, highC, Left - Acc, highC, right
        % Spd, lowC, Left - Spd, lowC, right - Spd, highC, Left - Spd, highC, right
    end
    
end

% Compute difference in cut-off RTs
diffT = allBestRTcut(subjects,useBins ~= 6) - repmat(allBestRTcut(subjects,useBins == 6),1,length(useBins)-1);
disp(['Mean signed difference in cut-off RT: '...
    num2str(nanmean(nanmean(diffT,2),1))...
    '+-'...
    num2str(nanstd(nanmean(diffT,2),1))]);
% Mean signed difference in cut-off RT: 0.0075286+-0.067508


%% Repeating the search for time of peak performance in individual 
% conditions to shows that the results are robust to this change

clear bestRTcond
binning = 3;
percentiles = [linspace(0,100,useBins(binning)+1)];

meanRT = nan(2,17,length(percentiles),2,2);
meanCorr = nan(2,17,length(percentiles),2,2);

for subj = subjects
    for cc=1:length(sats)
        for l = 1:2
            for lr = 1:2
                condsBehCA{cc,l,subj,1,lr} = find((indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                    indicators.respLR{subj}==indicators.LR{subj} &...
                    indicators.ContrLevels{subj} == l &...
                    indicators.LR{subj} == lr &...
                    indicators.RT{subj}>=minRT);
                condsBehCA{cc,l,subj,2,lr} = find((indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                    indicators.respLR{subj}~=indicators.LR{subj} &...
                    indicators.ContrLevels{subj} == l &...
                    indicators.LR{subj} == lr &...
                    indicators.RT{subj}>=minRT);
                
                thisTrials = [condsBehCA{cc,l,subj,1,lr} condsBehCA{cc,l,subj,2,lr}]; 
                numCorr = length(condsBehCA{cc,l,subj,1,lr});
                thisRT = indicators.RT{subj}(thisTrials);
                pc = [prctile(thisRT,percentiles)]; clear ts
                for p = 1:length(pc)-1
                    ts{p} = find(thisRT>=pc(p) & thisRT<=pc(p+1));
                    meanRT(lr,subj,p,cc,l) = nanmean(thisRT(ts{p})); 
                    meanCorr(lr,subj,p,cc,l) = length(find(ts{p}<=numCorr))./length(ts{p});
                end
            end
            
            thisCorr = squeeze(nanmean(nanmean(meanCorr(:,subj,:,cc,l),2),1));
            mpRT = find(thisCorr == max(thisCorr));
            bestRTcond(subj,cc,l) = mpRT(end);
            bestRTcut_cond(subj,cc,l) = squeeze(nanmean(meanRT(:,subj,bestRTcond(subj,cc,l),cc,l),1));
            
        end
    end
end

% Compute difference in cut-off RTs
diffT = bestRTcut_cond - repmat(allBestRTcut(:,binning),1,2,2);
disp(['Mean signed difference in cut-off RT: '...
    num2str(nanmean(nanmean(nanmean((diffT(subjects,:,:)),3),2)))...
    '+-'...
    num2str(nanstd(nanmean(nanmean((diffT(subjects,:,:)),3),2)))]);
% Mean signed difference in cut-off RT: 0.082864+-0.092802

disp('-------------------------------------------------')
disp('3-Way repeated-measures ANOVA on descencing slope of CAF past peak performance')
disp(['Cut-off RT based on 6 RT bins in indivudual conditions'])
cutRTs = bestRTcut_cond;
% CAF_slope_stats splits trials with RT > cutRTs into 8 RT bins, fits
% slopes to the data points within experimental conditions, computes
% 3-Way ANOVA on the slopes and prints the results
CAF_slope_stats(indicators, cutRTs, subjects);


% Speed/Accuracy: F(1,15)=0.59142, p=0.45381
% Contrast: F(1,15)=0.0387, p=0.84668
% Left/Right Target: F(1,15)=0.072337, p=0.79163
% Speed/Accuracy x Contrast: F(1,15)=0.01541, p=0.90286
% Speed/Accuracy x Left/Right Target: F(1,15)=0.068514, p=0.79707
% Contrast x Left/Right Target: F(1,15)=2.9326, p=0.1074
% Speed/Accuracy x Contrast x Left/Right Target: F(1,15)=0.14076, p=0.71278

