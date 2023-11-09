% Behaviour_stats

%% T-test and ANOVAs between speed and accuracy for RT and accuracy
s=0; clear meanRT meanAcc

for subj = subjects
    s=s+1; col = 0;
    
    for cc = 1:2    % 1 = Accuracy, 2 = Speed 
        
        for l = 1:2 % 1 = low Contrast, 2 = high Contrast
            
            col=col+1;
            
            meanRT(s,col) = nanmean(indicators.RT{subj}(find((indicators.cond{subj} == sats{cc}(1) | indicators.cond{subj} == sats{cc}(2)) &...
                indicators.ContrLevels{subj} == l)));
            
            allTr = length(find((indicators.cond{subj} == sats{cc}(1) | indicators.cond{subj} == sats{cc}(2)) &...
                indicators.ContrLevels{subj} == l));
            corrTr = length(find((indicators.cond{subj} == sats{cc}(1) | indicators.cond{subj} == sats{cc}(2)) &...
                indicators.ContrLevels{subj} == l &...
                indicators.LR{subj}==indicators.respLR{subj}));
            meanAcc(s,col) = corrTr/allTr;
        end
    end
end

% RT stats
output = teg_repeated_measures_ANOVA(meanRT, [2,2], {'Speed/Accuracy', 'Contrast'});
% allpsh = 0;
% [pWilk,pMauchly] = norm_sphere_test(meanRT);
% allpsh = allpsh+1; allSpheres.name{allpsh}='mean RT ANOVA'; allSpheres.pWilk{allpsh}=pWilk;
% allSpheres.pMauchly{allpsh}=pMauchly; save([resultsFolder '/allSpheres_behav'],'allSpheres')

disp('*************************************')
disp('RT')
for f = 1:size(output.R,1)
    disp([output.labels{f} ': F(1,15) = ' num2str(output.R(f,1)) ', p = ' num2str(output.R(f,4))])
end
disp('*************************************')

% Accuracy stats
output = teg_repeated_measures_ANOVA(meanAcc, [2,2], {'Speed/Accuracy', 'Contrast'});
% [pWilk,pMauchly] = norm_sphere_test(meanAcc);
% allpsh = allpsh+1; allSpheres.name{allpsh}='mean Accuracy ANOVA'; allSpheres.pWilk{allpsh}=pWilk;
% allSpheres.pMauchly{allpsh}=pMauchly; save([resultsFolder '/allSpheres_behav'],'allSpheres')

disp('*************************************')
disp('Response accuracy')
for f = 1:size(output.R,1)
    disp([output.labels{f} ': F(1,15) = ' num2str(output.R(f,1)) ', p = ' num2str(output.R(f,4))])
end
disp('*************************************')
% *************************************
% RT
% Speed/Accuracy: F(1,15) = 46.6268, p = 5.7114e-06
% Contrast: F(1,15) = 86.6709, p = 1.2697e-07
% Speed/Accuracy x Contrast: F(1,15) = 23.465, p = 0.00021455
% *************************************
% Response accuracy
% Speed/Accuracy: F(1,15) = 23.1752, p = 0.00022753
% Contrast: F(1,15) = 106.8142, p = 3.2331e-08
% Speed/Accuracy x Contrast: F(1,15) = 3.8607, p = 0.06823
% *************************************
%% Test for differences between points earned for speed and accuracy conditions

ordr = [1 4 3 2];
s=0;
for subj = subjects
    s=s+1; col=0;
    for col = 1:4
        trials = find(indicators.cond{subj} == ordr(col));
        meanRewANOVA(s,col) = nanmean(indicators.Rewards{subj}(1,trials));
    end
end

output = teg_repeated_measures_ANOVA(meanRewANOVA, [2,2], {'Speed/Accuracy', 'Deadline/Slope'});
% [pWilk,pMauchly] = norm_sphere_test(meanRewANOVA);
% allpsh = allpsh+1; allSpheres.name{allpsh}='mean Reward ANOVA'; allSpheres.pWilk{allpsh}=pWilk;
% allSpheres.pMauchly{allpsh}=pMauchly; save([resultsFolder '/allSpheres_behav'],'allSpheres')

disp('*************************************')
disp('Reward')
for f = 1:size(output.R,1)
    disp([output.labels{f} ': F(1,15) = ' num2str(output.R(f,1)) ', p = ' num2str(output.R(f,4))])
end
disp('*************************************')

% *************************************
% Reward
% Speed/Accuracy: F(1,15) = 68.8487, p = 5.4765e-07
% Deadline/Slope: F(1,15) = 4.4866, p = 0.051271
% Speed/Accuracy x Deadline/Slope: F(1,15) = 64.4537, p = 8.2442e-07
% *************************************


savepath = fullfile(dataFolder, 'ANOVA data for SPSS');
if exist(savepath, 'dir') ~= 7
    mkdir(savepath);
end
save(fullfile(savepath, 'behavior_anovas'),...
    'meanRT','meanAcc','meanRewANOVA')
%% Miss rates across subjects for Speed/Accuracy/Total
% Trials during which the subject did not resopnd before the offset of the
% stimulus were classified as 'missed'

missRate = nan(17,2,2); missRateAll = nan(1,17);

for cc = 1:2 % 1 = Accuracy, 2 = Speed
    for l = 1:2 % 1 = low Contrast, 2 = high Contrast
        for subj = subjects
            
            trials = find((indicators.cond{subj} == sats{cc}(1) | indicators.cond{subj} == sats{cc}(2)) &...
                indicators.ContrLevels{subj} == l);
            
            trialsNoR = find((indicators.cond{subj} == sats{cc}(1) | indicators.cond{subj} == sats{cc}(2)) &...
                indicators.ContrLevels{subj} == l &...
                indicators.respTypes{subj} == 5); % 'No response' trials were marked as category 5 during the experiment
            
            missRate(subj,cc,l) = length(trialsNoR) / length(trials);
            
        end
        disp([AccSpd{cc} ', ' LoHi{l} ': miss rate = ' num2str(nanmean(missRate(subjects,cc,l))) ' + - ' num2str(nanstd(missRate(subjects,cc,l),[],1))])
    end
end

disp(['Overall miss rate = ' num2str(nanmean(nanmean(nanmean(missRate(subjects,:,:))))) ' + - ' num2str(nanstd(nanmean(nanmean(missRate(subjects,:,:),2),3),[],1))])
% Accuracy, Low Contrast: miss rate = 0.0034999 + - 0.0080657
% Accuracy, High Contrast: miss rate = 0.00026824 + - 0.001073
% Speed, Low Contrast: miss rate = 0.00080244 + - 0.0023307
% Speed, High Contrast: miss rate = 0.00079903 + - 0.0017179
% Overall miss rate = 0.0013424 + - 0.0025906

%% Effect of Speed/Accuracy on decrease in response accuracy over RT
% Measure slope of descending conditional accuracy past the point of peak
% performance to determine whether the drop in performance is greater under
% speed pressure 
% (This is a script, not a function)

CAF_tail_stats

%% Plots and statistical analyses related to Deadline and Slope conditions 
% within the Speed and Accuracy Regime can be found in SuppFigure1_plot.m 





