% SuppFigure5_plot



disp('Loading STFT SL...'); load(fullfile(figData, ['CSDstft_' conditions{1} '.mat']));

% Mu/Beta amplitude just before evidence onset (window: -300 to 0ms
% before evidence onset) split into RT-bins plotted over RT

% Compute traces for all four conditions

clear condsMuPooledBinnedM1
for subj = subjects
    for c=1:max(indicators.cond{subj})
        for l = 1:2
            for d = 1:3
                for lr = 1:2
                    clear trlDelay rtlimDelay
                    trlDelay = find(indicators.onsedelay{subj} == d  &...
                        indicators.cond{subj}==c &...
                        indicators.ContrLevels{subj} == l &...
                        indicators.LR{subj} == lr &...
                        validrlockS{subj} &...
                        goodTrialsComb{subj});
                    
                    rtlimDelay = prctile(indicators.RT{subj}(trlDelay),[linspace(0,100,CPPnumbins+1)]);
                    for b = 1:CPPnumbins
                        condsMuPooledBinnedM1{c,l,d,subj,lr,b} = find(indicators.RT{subj}>=rtlimDelay(b) &...
                            indicators.RT{subj}<rtlimDelay(b+1) &...
                            indicators.onsedelay{subj} == d &...
                            indicators.cond{subj}==c &...
                            indicators.ContrLevels{subj} == l &...
                            indicators.respLR{subj} == lr &...
                            validrlockS{subj} &...
                            goodTrialsComb{subj});
                    end
                    
                    clear trlDelay rtlimDelay
                    trlDelay = find(indicators.onsedelay{subj} == d &...
                            indicators.cond{subj}==c &...
                            indicators.ContrLevels{subj} == l &...
                            indicators.respLR{subj} == lr &...
                            indicators.respLR{subj} == indicators.LR{subj} &...
                            validrlockS{subj} &...
                            goodTrialsComb{subj});
                        
                    rtlimDelay = prctile(indicators.RT{subj}(trlDelay),[linspace(0,100,CPPnumbins+1)]);
                    for b = 1:CPPnumbins
                        condsMuCorrectBinned{c,l,d,subj,lr,b} = find(indicators.RT{subj}>=rtlimDelay(b) &...
                            indicators.RT{subj}<rtlimDelay(b+1) &...
                            indicators.onsedelay{subj} == d &...
                            indicators.cond{subj}==c &...
                            indicators.ContrLevels{subj} == l &...
                            indicators.respLR{subj} == lr &...
                            indicators.respLR{subj} == indicators.LR{subj} &...
                            validrlockS{subj} &...
                            goodTrialsComb{subj});
                    end
                    
                end
            end
        end
    end
end

meanMu_BL = NaN*ones(CPPnumbins,2,2,17,3,2,2); 
meanMu_Exc = NaN*ones(CPPnumbins,2,2,17,3,2); 
meanRTCtr = NaN*ones(CPPnumbins,2,2,17,3,2);
meanRT_Exc = NaN*ones(CPPnumbins,2,2,17,3,2);
for cc=1:2
    for l = 1:2
        for b=1:CPPnumbins
            for subj = subjects
                for d = 1:3
                    for lr = doLR
                        trials = [condsMuPooledBinnedM1{sats{cc}(1),l,d,subj,lr,b} condsMuPooledBinnedM1{sats{cc}(2),l,d,subj,lr,b}];
                        meanRTCtr(b,cc,l,subj,d,lr) = nanmean(indicators.RT{subj}(1,trials)*1000,2)';
                        
                        meanMu_BL(b,cc,l,subj,d,lr,1) = squeeze(nanmean(nanmean(CSDstft{subj}(mueChans(3-lr),Frange,find(T == MBrBin),trials),4),2));
                        meanMu_BL(b,cc,l,subj,d,lr,2) = squeeze(nanmean(nanmean(CSDstft{subj}(mueChans(lr),Frange,find(T == MBrBin),trials),4),2));
                        
                        trials = [condsMuCorrectBinned{sats{cc}(1),l,d,subj,lr,b} condsMuCorrectBinned{sats{cc}(2),l,d,subj,lr,b}];
                        meanMu_Exc(b,cc,l,subj,d,lr) = squeeze(nanmean(nanmean(CSDstftr{subj}(mueChans(lr),Frange,find(Tr == MBrBin),trials),4),2)) -...
                            squeeze(nanmean(nanmean(CSDstft{subj}(mueChans(lr),Frange,find(T == MBrBin),trials),4),2));
                        meanRT_Exc(b,cc,l,subj,d,lr) = nanmean(indicators.RT{subj}(1,trials)*1000,2)';
                    end
                end
            end
        end
    end
end


save(fullfile(figData, 'MB_BL_over_RT'), 'meanRTCtr', 'meanMu_BL', 'meanMu_Exc', 'meanMu_Exc')

