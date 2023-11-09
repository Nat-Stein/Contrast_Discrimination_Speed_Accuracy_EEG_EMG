% SuppFigure3_prepstats


%% -----------------------------------------------------------------
%% 2. CPP slope stats using ANOVAs
% Compute rate of rise in CPP amplitude (with Auditory Evoked Potential)
% Build-up of CPP scales with Contrast and is increased under speed pressure

% load(fullfile(figData, ['CSDerpr_' conditions{1} '.mat']),'CSDerpr')

N = 4; [B,A] = butter(N,8*2/fs);
CPPSlope_AEP = NaN*ones(max(subjects),2,2,3,2);

tslope = find(tr >= -300 & tr <= -50); % This is the final CPP slope time frame (18-04-02)

for subj = subjects
    CPPslopeSaep{subj} = nan(1,size(CSDerpr{subj},3));
    cfSlope = nan(1,size(CSDerpr{subj},3));
    
    for tt = 1 : length(indicators.RT{subj})
        
        % Fit line to centro-parietal positivity in time window between 300
        % and 50 ms before the button click
        
        inDat = nanmean(CSDerpr{subj}(chCPP,:,tt),3);
        
        if sum(isnan(inDat)) ~= length(inDat)
            this = filtfilt(B, A, inDat);
            [p,S] = polyfit(tr(tslope)-tr(tslope(1)),this(:,tslope),1);
            cfSlope(1,tt) = p(1);
        else
            cfSlope(1,tt) = nan;
        end
    end
    CPPslopeSaep{subj} = cfSlope;
    
    for cc = 1:2
        for l = 1:2
            for d = 1:3
                for lr = 1:2
                    
                    trials = find(indicators.onsedelay{subj} == d  &...
                        (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                        indicators.ContrLevels{subj} == l &...
                        indicators.LR{subj} == lr &...
                        validrlockS{subj} &...
                        goodTrialsComb{subj});
                    
                    
                    CPPSlope_AEP(subj,cc,l,d,lr) = nanmean(CPPslopeSaep{subj}(1,trials));
                    
                end
            end
        end
    end
end

save(fullfile(figData, 'CPPSlope_AEP'),'CPPSlope_AEP')
save(fullfile(figData, 'CPPslopeSaep'),'CPPslopeSaep','tslope')


































