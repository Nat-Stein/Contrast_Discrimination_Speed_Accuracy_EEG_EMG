%% Figure4_prepstats


%% 1. Build data table for Mu/Beta amplitude at response

load(fullfile(figData, ['CSDstftr_' conditions{1} '.mat']))

MBrBin = -fftlen;
MBatRT = []; RT = []; SAT = []; Contrast = []; Subject = []; IpsCon = []; LvsR = [];

satConds = [1 2 2 1];
for subj = subjects
    clear allTrials thisMB allMB allRT
    allMB = nan(2,length(indicators.RT{subj})); allRT = nan(1,length(indicators.RT{subj}));
    
    allTrials = find((goodTrialsComb{subj} == 1)' &...
        ~isnan(squeeze(nanmean(CSDstftr{subj}(mueChans(lr),Frange,find(Tr == MBrBin),:),2))) &...
        ~isnan(indicators.RT{subj})');
    
    allMB(1,allTrials) = zscore(nanmean(CSDstftr{subj}(mueChans(1),Frange,find(Tr == MBrBin),allTrials),2));
    allMB(2,allTrials) = zscore(nanmean(CSDstftr{subj}(mueChans(2),Frange,find(Tr == MBrBin),allTrials),2));
    allRT(1,allTrials) = zscore(indicators.RT{subj}(allTrials));
    
    for lr = 1:2
        
        trials = find(goodTrialsComb{subj} == 1 &...
            squeeze(~isnan(nanmean(CSDstftr{subj}(mueChans(lr),Frange,find(Tr==MBrBin),:),2)))' &...
            indicators.respLR{subj} == lr);
        
        for ic = 1:2 % 1 = contralateral to executed response, 
            % 2 = ipsilateral to executed response
            
            if ic == 1
                thisMB = squeeze(allMB(3-lr,trials)); % Contralalteral
            else
                thisMB = squeeze(allMB(lr,trials));
            end
            
            RTS = allRT(1,trials);
            
            MBatRT = [MBatRT thisMB];
            RT = [RT RTS];
            SAT = [SAT satConds(indicators.cond{subj}(trials))];
            Contrast = [Contrast indicators.ContrLevels{subj}(trials)];
            Subject = [Subject subj*ones(1,length(trials))];
            
            IpsCon = [IpsCon ic*ones(1,length(trials))];
            LvsR = [LvsR lr*ones(1,length(trials))];
        end
    end
end

save(fullfile(figData, 'MB_rtLME.mat'), 'MBatRT', 'Subject', 'RT', 'SAT', 'Contrast', 'LvsR', 'IpsCon')

%% ------------------------------------------------------------------
%% 2. ANOVA for MuBeta amplitude at response: executed vs withheld response

% Run if CSDstftr is not already loaded
load(fullfile(figData, ['CSDstftr_' conditions{1} '.mat']))

MBrBin = -fftlen;

corrOnly = 0;

meanMBresp = nan(17,2,2,2); meanMBrespNoZ = nan(17,2,2,2);
satConds = [1 2 2 1];
for subj = subjects
    clear allTrials thisMB allMB allRT
    allMB = nan(2,length(indicators.RT{subj})); allRT = nan(1,length(indicators.RT{subj}));
    
    allTrials = find((goodTrialsComb{subj}==1)' &...
        ~isnan(squeeze(nanmean(CSDstftr{subj}(mueChans(lr),Frange,find(Tr==MBrBin),:),2))) &...
        ~isnan(indicators.RT{subj})');
    
    allMB(1,allTrials) = zscore(nanmean(CSDstftr{subj}(mueChans(1),Frange,find(Tr==MBrBin),allTrials),2));
    allMB(2,allTrials) = zscore(nanmean(CSDstftr{subj}(mueChans(2),Frange,find(Tr==MBrBin),allTrials),2));
    
    for cc = 1:2
        for l = 1:2
            for ic = 1:2
                thisMB = []; thismb = [];
                for lr = 1:2
                    trials = find(goodTrialsComb{subj}==1 &...
                        squeeze(~isnan(nanmean(CSDstftr{subj}(mueChans(lr),Frange,find(Tr==MBrBin),:),2)))' &...
                        indicators.respLR{subj} == lr &...
                        (indicators.cond{subj} == sats{cc}(1) | indicators.cond{subj} == sats{cc}(2)) &...
                        indicators.ContrLevels{subj} == l);
                    
                    if ic == 1
                        thisMB = [thisMB squeeze(allMB(3-lr,trials))];
                    else
                        thisMB = [thisMB squeeze(allMB(lr,trials))];
                    end
                end
                meanMBresp(subj,cc,l,ic) = nanmean(thisMB);
            end
        end
    end
end

save(fullfile(figData, 'MB_rtLME_executed_vs_withheld.mat'), 'meanMBresp')


% ------------------------------------------------------------------
%% 3. CPP amplitude at baseline (before evidence onset)

load(fullfile(figData, ['CSDerp_' conditions{2} '.mat']))
cue = load(fullfile(figData, ['goodTrials__' conditions{2}])); goodTrialsCombCue = cue.goodTrialsComb; clear cue
load(fullfile(figData, ['indicators_' conditions{2}]));

corrOnly = 0;
t_BL = [-50 -2];
CPPpure = []; CPPride = []; RT = []; SAT = []; Contrast = []; Subject = []; LvsR = [];

satConds = [1 2 2 1];
for subj = subjects
    
        trials = find(goodTrialsCombCue{subj}==1 &...
            squeeze(~isnan(nanmean(CSDerp{subj}(chCPP,find(t>=t_BL(1) & t<=t_BL(2)),:),2)))');
    
    CPPrideS = zscore(squeeze(nanmean(CSDerp{subj}(chCPP,find(t>=t_BL(1) & t<=t_BL(2)),trials),2)))';
    RTS = zscore(indicatorsCue.RT{subj}(trials));
    
    CPPride = [CPPride CPPrideS];
    RT = [RT RTS];
    SAT = [SAT satConds(indicatorsCue.cond{subj}(trials))];
    Contrast = [Contrast indicatorsCue.ContrLevels{subj}(trials)];
    Subject = [Subject subj*ones(1,length(trials))];
    LvsR = [LvsR indicatorsCue.LR{subj}(trials)];
end

save(fullfile(figData, 'CPP_blLME.mat'),'CPPride','Subject','RT','SAT','Contrast','LvsR','MBrBin','t_BL')


% ------------------------------------------------------------------
%% 4. CPP amplitude at response

for calcWithAEP = 0 : 1
    
    if calcWithAEP == 0  % Calculate for RIDE-corrected ERPs with AEP removed
        load(fullfile(figData, 'CSDerprRIDE'))
    else % Calculate for original ERPs with AEP
        load(fullfile(figData, ['CSDerpr_' conditions{1} '.mat']),'CSDerpr')
    end
    
    CPPride = []; CPPaep = []; RT = []; SAT = []; Contrast = []; Subject = []; LvsR = [];
    
    t_RT = [-130 -70];
    satConds = [1 2 2 1];
    for subj = subjects
        
        if calcWithAEP == 1
            trials = find(goodTrialsComb{subj}==1 &...
                squeeze(~isnan(nanmean(CSDerpr{subj}(chCPP,find(tr>=t_RT(1) & tr<t_RT(2)),:),2)))');
            CPPsubj = zscore(squeeze(nanmean(CSDerpr{subj}(chCPP,find(tr>=t_RT(1) & tr<t_RT(2)),trials),2)))';
            CPPaep = [CPPaep CPPsubj];
        else
            trials = find(goodTrialsComb{subj}==1 &...
                squeeze(~isnan(nanmean(CSDerprRIDE{subj}(chCPP,find(tr>=t_RT(1) & tr<t_RT(2)),:),2)))');
            CPPsubj = zscore(squeeze(nanmean(CSDerprRIDE{subj}(chCPP,find(tr>=t_RT(1) & tr<t_RT(2)),trials),2)))';
            CPPride = [CPPride CPPsubj];
        end
        
        RT = [RT zscore(indicators.RT{subj}(trials))];
        SAT = [SAT satConds(indicators.cond{subj}(trials))];
        Contrast = [Contrast indicators.ContrLevels{subj}(trials)];
        Subject = [Subject subj*ones(1,length(trials))];
        LvsR = [LvsR indicators.LR{subj}(trials)];
    end
    
    
    if calcWithAEP == 0
        
        datset=dataset({CPPride','CPPride'},...
            {Subject','Subject'},...
            {RT','RT'},...
            {SAT','SAT'},...
            {Contrast','Contrast'},...
            {LvsR','LvsR'});
        
        save(fullfile(figData, 'CPP_at_RT_stats'),'datset','fixedE','randomE','t_RT')
        
    else
        
        datset=dataset({CPPaep','CPPaep'},...
            {Subject','Subject'},...
            {RT','RT'},...
            {SAT','SAT'},...
            {Contrast','Contrast'},...
            {LvsR','LvsR'});
        
        save(fullfile(figData, 'CPP_at_RT_stats_AEP'),'datset','fixedE','randomE','t_RT')
        
    end
    
end


%% 5. CPP slope stats using ANOVAs
% Build-up of CPP scales with Contrast and is increased under speed pressure

load(fullfile(figData, 'CSDerprRIDE'))

N = 4; [B,A] = butter(N,8*2/fs);
CPPSlopeRIDE = NaN*ones(max(subjects),2,2,3,2);

tslope = find(tr >= -300 & tr <= -50); % This is the final CPP slope time frame (18-04-02)

for subj = subjects
    CPPslopeS{subj} = nan(1,size(CSDerprRIDE{subj},3));
    cfSlopeRIDE = nan(1,size(CSDerprRIDE{subj},3));
    
    for tt = 1:length(indicators.RT{subj})
        
        % Fit line to centro-parietal positivity in time window between 300
        % and 50 ms before the button click
        e = nanmean(CSDerprRIDE{subj}(chCPP,:,tt),3);
        if sum(isnan(e)) < length(e)
            this = filtfilt(B,A,nanmean(CSDerprRIDE{subj}(chCPP,:,tt),3));
            [p,S] = polyfit(tr(tslope)-tr(tslope(1)),this(tslope),1);
            
            cfSlopeRIDE(1,tt) = p(1);
        else
            cfSlopeRIDE(1,tt) = nan;
        end
    end
    CPPslopeS{subj} = cfSlopeRIDE;
    
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

                    CPPSlopeRIDE(subj,cc,l,d,lr) = nanmean(cfSlopeRIDE(1,trials));
                    
                end
            end
        end
    end
end

save(fullfile(figData, 'CPPSlopeRIDE'), 'CPPSlopeRIDE')
save(fullfile(figData, 'CPPslopeS'), 'CPPslopeS','tslope')


% ------------------------------------------------------------------
%% 6. Mu/Beta slope stats using ANOVAs

% If not loaded already
% disp('Loading STFT RL...'); load(fullfile(figData, ['CSDstftr_' conditions{1}]));

% Prepare data
tMB = find(Tr >= -350 & Tr <= -150);

MBSlope = NaN*ones(max(subjects),2,2,3,2,2);

for subj = subjects
    MBslopeS{subj} = nan(size(CSDstftr{subj}));
    mbSlope = NaN*ones(2,size(CSDstftr{subj},4));
    
    for tt = 1:size(mbSlope,2)
        mbThis = squeeze(nanmean(nanmean(CSDstftr{subj}(mueChans,Frange,:,tt),4),2));
        [p,S] = polyfit([1:length(tMB)],squeeze(mbThis(1,tMB)),1); mbSlope(1,tt) = p(1);
        [p,S] = polyfit([1:length(tMB)],squeeze(mbThis(2,tMB)),1); mbSlope(2,tt) = p(1);
    end
    MBslopeS{subj} = mbSlope;
    
    for cc = 1:2
        for l = 1:2
            for d = 1:3
                for lr = 1:2
                    
                    trials = find(indicators.LR{subj} == indicators.respLR{subj} & indicators.onsedelay{subj} == d  &...
                        (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                        indicators.ContrLevels{subj} == l &...
                        indicators.respLR{subj} == lr & validrlockS{subj} & goodTrialsComb{subj});
                    
                    MBSlope(subj,cc,l,d,lr,1) = nanmean(MBslopeS{subj}((3-lr),trials));
                    MBSlope(subj,cc,l,d,lr,2) = nanmean(MBslopeS{subj}(lr,trials));
                end
            end
        end
    end
end

save(fullfile(figData, 'MBslopeS'),'MBslopeS','mbSlope_anova');






