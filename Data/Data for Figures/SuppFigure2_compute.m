% SuppFigure2_compute.m
% Calculate SSVEP frequencies over time as well as other frequencies for
% comparison

% load(fullfile(dataFolder, conditions{1}, 'SSVEP_freqSpectrum')) 
load(fullfile(figData, 'SSVEP_freqSpectrum')) 
% Computed in Figure2_compute


meanSL_Freq = NaN*ones(2,2,2,length(Fssvep),length(T));
meanSL_FreqSubj = NaN*ones(2,2,2,length(Fssvep),length(T),max(subjects));
RTmeanS = nan(2,2,2,17);

for cc = 1:2
    for l = 1:2
        for lr = doLR
            thisSSVEP = NaN*ones(3,17,length(Fssvep),length(T));
            RTmean = nan(3,17);
            for subj = subjects
                for d = 1:3
                    
                    trialsSSVEP = find(indicators.onsedelay{subj} == d  &...
                        (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                        indicators.ContrLevels{subj} == l &...
                        indicators.LR{subj} == lr &...
                        validrlockS{subj} &...
                        goodTrialsComb{subj});
                    
                    % SSVEP
                    thisSSVEP(d,subj,:,:) = squeeze(nanmean(allFreq{subj}(:,:,trialsSSVEP),3));
                    
                    % Mean RTs
                    RTmean(d,subj) = nanmean(indicators.RT{subj}(trialsSSVEP));
                end
                
                % SSVEP
                meanSL_FreqSubj(cc,l,lr,:,:,subj) = squeeze(nanmean(thisSSVEP(d,subj,:,:),1));
                
                % Mean RTs
                RTmeanS(cc,l,lr,subj) = nanmean(RTmean(:,subj));
            end
        end
    end
end

save(fullfile(dataFolder, 'SSVEPindFreq'),'meanSL_FreqSubj','RTmeanS')


%% Calculate stats for whole time window that's significant for the 
% differntial SSVEP

twin = [350 550]; 

take_freqs = [20 22.5 25];
for freq = 1:length(take_freqs)
    thisF = find(Fssvep == take_freqs(freq));
    clear freq_anova
    con=0;
    for cc=1:2
        for l = 1:2
            for lr = doLR
                s=0; con = con+1;
                for subj = subjects
                    s=s+1;
                    freq_anova(s,con) = squeeze(nanmean(meanSL_FreqSubj(cc,l,lr,thisF,find(T>=twin(1) & T<=twin(2)),subj,1),5));
                end
            end
        end
    end
    output = teg_repeated_measures_ANOVA(freq_anova, [2,2,2], {'SAT','Contrast','LR'});
    
    windowPs(freq,:) = output.R(:,4);
    windowFs(freq,:) = output.R(:,1);
    disp('-----------------')
    disp([num2str(take_freqs(freq)) 'Hz'])
    for f = 1:size(output.R,1)
        disp([output.labels{f} ': F(1,15) = ' num2str(output.R(f,1)) ', p = ' num2str(output.R(f,4))])
    end
    
    % Save data for Bayes Factor analysis in JASP
    regression_table_singleFreq = [[1:size(freq_anova,2)]; freq_anova];
    csvwrite(fullfile(dataFolder, ['regression_table_singleFreq' num2str(take_freqs(freq)) '.csv']),regression_table_singleFreq)
    % Format:
    % Acc, lowC, Left - Acc, lowC, right - Acc, highC, Left - Acc, highC, right
    % Spd, lowC, Left - Spd, lowC, right - Spd, highC, Left - Spd, highC, right
    
end

% 20Hz: SAT: F(1,15) = 0.36178, p = 0.5565
% 22.5Hz SAT: F(1,15) = 0.0058297, p = 0.94015
% 25Hz SAT: F(1,15) = 2.4941, p = 0.13512
% -----------------
% 20Hz
% SAT: F(1,15) = 0.36178, p = 0.5565
% Contrast: F(1,15) = 24.2627, p = 0.00018295
% LR: F(1,15) = 8.9762, p = 0.0090455
% SAT x Contrast: F(1,15) = 2.3902, p = 0.14293
% SAT x LR: F(1,15) = 0.52957, p = 0.47799
% Contrast x LR: F(1,15) = 8.6487, p = 0.010119
% SAT x Contrast x LR: F(1,15) = 1.9378, p = 0.1842
% -----------------
% 22.5Hz
% SAT: F(1,15) = 0.0058297, p = 0.94015
% Contrast: F(1,15) = 0.086579, p = 0.7726
% LR: F(1,15) = 1.7539, p = 0.20521
% SAT x Contrast: F(1,15) = 1.727, p = 0.20854
% SAT x LR: F(1,15) = 0.77471, p = 0.39265
% Contrast x LR: F(1,15) = 4.0573, p = 0.062276
% SAT x Contrast x LR: F(1,15) = 1.8912, p = 0.18926
% -----------------
% 25Hz
% SAT: F(1,15) = 2.4941, p = 0.13512
% Contrast: F(1,15) = 1.5951, p = 0.22588
% LR: F(1,15) = 21.9351, p = 0.00029415
% SAT x Contrast: F(1,15) = 0.71755, p = 0.41026
% SAT x LR: F(1,15) = 1.9676, p = 0.18106
% Contrast x LR: F(1,15) = 23.5782, p = 0.00020971
% SAT x Contrast x LR: F(1,15) = 0.75392, p = 0.39892


%% Compute significance of main effect of Speed/Accuracy emphasis in 
% every time window

take_freqs = [20 22.5 25];
allPs = nan(length(take_freqs),length(T),7);
allFs = nan(length(take_freqs),length(T),7);
for freq = 1:length(take_freqs)
    thisF = find(Fssvep == take_freqs(freq));
    for tpt = 1:length(T)
        clear freq_anova
        con=0; 
        for cc=1:2
            for l = 1:2
                for lr = doLR
                    s=0; con = con+1;
                    for subj = subjects
                        s=s+1;
                        freq_anova(s,con) = squeeze(meanSL_FreqSubj(cc,l,lr,thisF,tpt,subj));
                    end
                end
            end
        end
        output = teg_repeated_measures_ANOVA(freq_anova, [2,2,2], {'SAT','Contrast','LR'});
        
        allPs(freq,tpt,:) = output.R(:,4);
        allFs(freq,tpt,:) = output.R(:,1);
    end
end

% save(fullfile(dataFolder, 'SSVEPindFreq_F'),'allFs','allPs')
save(fullfile(figData, 'SSVEPindFreq_F'),'allFs','allPs')






