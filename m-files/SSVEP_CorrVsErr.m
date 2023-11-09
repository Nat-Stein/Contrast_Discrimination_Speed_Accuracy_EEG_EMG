% SSVEP_CorrVsErr.m: Test SSVEP behavioral effect

% Load single-trial SSVEPs computed in Figure2_compute.m
% SSVEPsubtr is the SSVEP minus it's neighbours
% Structure SSVEPsubtr{subject}(frequency,time point,trial number)

load(fullfile(figData, 'SSVEPsubtr_only'),'SSVEPsubtr');
load(fullfile(figData, 'indicators_stimLock'))


%% Calculate difference between SSVEPs

% Calculate SSVEP power difference between increased and decreased 
% frequency per trial and time point
for subj = subjects
    meanCorrSSVEP{subj} = nan(size(SSVEPsubtr{subj},2),size(SSVEPsubtr{subj},3));
    for lr = 1:2
        for l = 1:2
            trialsSSVEP = find(indicators.LR{subj} == lr &...
                indicators.ContrLevels{subj} == l); 
            % indicators.LR{subj} specifies which frequency increased
            % indicators.ContrLevels{subj} gives contrast difference level
            
            clear thisSSVEP; thisSSVEP = squeeze(SSVEPsubtr{subj}(lr,:,trialsSSVEP) - SSVEPsubtr{subj}(3-lr,:,trialsSSVEP));
            meanCorrDiffSSVEP{subj}(:,trialsSSVEP) = thisSSVEP;
            diffSSVEP_lr{subj}(:,trialsSSVEP) = squeeze(SSVEPsubtr{subj}(1,:,trialsSSVEP) - SSVEPsubtr{subj}(2,:,trialsSSVEP));

        end
    end
end

meanResSSVEP_Diff_LR = nan(max(subjects),2,size(meanCorrSSVEP{subj},1),2,2);
for subj = subjects % Subjects
    for l = 1:2 % Contrast levels: 1 = low, 2 = high       
        for lr = 1:2 % left (1) vs right (2) contrast increase
            % Take mean of all CORRECT trials with same contrast & side
            % indicators.respLR{subj} gives the subjects response
            trials = find(indicators.LR{subj} == lr & indicators.LR{subj} == indicators.respLR{subj} &...
                indicators.ContrLevels{subj} == l);
            meanResSSVEP_Diff_LR(subj,1,:,l,lr) = squeeze(nanmean(meanCorrDiffSSVEP{subj}(:,trials),2));
            diffSSVEP_lr_LvsR(subj,1,:,l,lr) = squeeze(nanmean(diffSSVEP_lr{subj}(:,trials),2));
            
            % Take mean of all WRONG trials with same contrast & side
            trials = find(indicators.LR{subj} == lr & indicators.LR{subj} ~= indicators.respLR{subj} &...
                indicators.ContrLevels{subj} == l);
            meanResSSVEP_Diff_LR(subj,2,:,l,lr) = squeeze(nanmean(meanCorrDiffSSVEP{subj}(:,trials),2));   
            diffSSVEP_lr_LvsR(subj,2,:,l,lr) = squeeze(nanmean(diffSSVEP_lr{subj}(:,trials),2));   
        end
    end
end


%% Plot left-right over time taken from SSVEP (left) minus SSVEP(right)
% 3-Way ANOVA with main effects of Correct/Wrong, Left/Right trial,
% Speed/Accuracy

t_bin = [350 550];
% t_bin gives the range of time points which will be averaged to then
% calculate an ANOVA
% The time between 200 and 550ms after evidence-onset was chosen, because
% the effect of Speed/Accuracy Regime on SSVEP was significant for all time
% point in that window
clear ssvepDiffRes_anova thisPlot; 
con = 0;
ssvepDiff_lr_sat_anova = nan(16,8);
for l = 1 
    % We are only considering low Contrast trials, because most subjects have very few errors on high contrast trials
    for cw = 1:2 % correct (1) vs wrong (2)
        for lr = 1:2 % 1 = left trial, 2 = right trial
            for cc = 1:2
                s=0; con = con+1;
                for subj = subjects
                    s=s+1;
                    if cw == 1
                        % Correct response trials
                        trials = find(indicators.LR{subj} == lr &...
                            indicators.LR{subj} == indicators.respLR{subj} &...
                            indicators.ContrLevels{subj} == l &...
                            (indicators.cond{subj} == sats{cc}(1) | indicators.cond{subj} == sats{cc}(2)));
                    else
                        % Incorrect response trials
                        trials = find(indicators.LR{subj} == lr &...
                            indicators.LR{subj} ~= indicators.respLR{subj} &...
                            indicators.ContrLevels{subj} == l &...
                            (indicators.cond{subj} == sats{cc}(1) | indicators.cond{subj} == sats{cc}(2)));
                    end
                    
                    % Prepare matrix for ANOVA
                    ssvepDiff_lr_sat_anova(s,con) = squeeze(nanmean(nanmean(diffSSVEP_lr{subj}(find(T>= t_bin(1) & T<= t_bin(2)),trials))));
                end
            end
        end
    end
end

%% Compute ANOVA on mean SSVEP per Correct/Wrong and Left/Right

output = teg_repeated_measures_ANOVA(ssvepDiff_lr_sat_anova, [2,2,2], {'Correct/Wrong','LvsR','SAT'});
for f = 1:size(output.R,1)
    disp([output.labels{f} ': F(1,15) = ' num2str(output.R(f,1)) ', p = ' num2str(output.R(f,4))])
end

% Correct/Wrong: F(1,15) = 0.38334, p = 0.54511
% LvsR: F(1,15) = 21.4874, p = 0.00032345
% SAT: F(1,15) = 0.029145, p = 0.86673
% Correct/Wrong x LvsR: F(1,15) = 7.7059, p = 0.014128
% Correct/Wrong x SAT: F(1,15) = 0.036238, p = 0.85158
% LvsR x SAT: F(1,15) = 4.6193, p = 0.048342
% Correct/Wrong x LvsR x SAT: F(1,15) = 0.072853, p = 0.79091





