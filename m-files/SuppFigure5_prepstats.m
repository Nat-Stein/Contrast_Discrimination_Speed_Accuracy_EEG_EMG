% SuppFigure5_prepstats




load(fullfile(figData, ['CSDstft_' conditions{1} '.mat']))

MBrBin = -fftlen;
MBatBL = []; RT = []; SAT = []; Contrast = []; Subject = []; IpsCon = []; LvsR = []; MBdiff = []; isCorr = [];

corrOnly = 0; addContrast = 0;
satConds = [1 2 2 1];
for subj = subjects
    clear allTrials thisMB allMB allRT
    allMB = nan(2,length(indicators.RT{subj})); allRT = nan(1,length(indicators.RT{subj}));
    allTrials = find(~isnan(squeeze(nanmean(CSDstft{subj}(mueChans(lr),Frange,find(T==MBrBin),:),2))) & ~isnan(indicators.RT{subj})');
    allMB(1,allTrials) = zscore(nanmean(CSDstft{subj}(mueChans(1),Frange,find(T==MBrBin),allTrials),2));
    allMB(2,allTrials) = zscore(nanmean(CSDstft{subj}(mueChans(2),Frange,find(T==MBrBin),allTrials),2));
    allRT(1,allTrials) = zscore(indicators.RT{subj}(allTrials));
    
    for lr = 1:2
        
        if corrOnly == 1
            trials = find(goodTrialsComb{subj}==1 &...
                squeeze(~isnan(nanmean(CSDstftr{subj}(mueChans(lr),Frange,find(T==MBrBin),:),2)))' &...
                indicators.respLR{subj}==indicators.LR{subj} & indicators.respLR{subj} == lr);
        else
            trials = find(goodTrialsComb{subj}==1 &...
                squeeze(~isnan(nanmean(CSDstftr{subj}(mueChans(lr),Frange,find(T==MBrBin),:),2)))' &...
                indicators.respLR{subj} == lr);
        end
        
        for ic = 1:2
            
            if ic == 1
                thisMB = squeeze(allMB(3-lr,trials));
                thisDiff = squeeze(allMB(3-lr,trials)) - squeeze(allMB(lr,trials));
            else
                thisMB = squeeze(allMB(lr,trials));
                thisDiff = nan(1,length(trials));
            end
            RTS = allRT(1,trials);
            
            MBdiff = [MBdiff thisDiff];
            MBatBL = [MBatBL thisMB];
            RT = [RT RTS];
            SAT = [SAT satConds(indicators.cond{subj}(trials))];
            Subject = [Subject subj*ones(1,length(trials))];
            Contrast = [Contrast indicators.ContrLevels{subj}(trials)];
            
            isCorr = [isCorr indicators.respLR{subj}(trials) == indicators.LR{subj}(trials)];
            IpsCon = [IpsCon ic*ones(1,length(trials))];
            LvsR = [LvsR lr*ones(1,length(trials))];
        end
    end
end



%% Fitlme for MB at baseline
exec = {'executed','withheld'};

% Mu/Beta activity is measured in the 300-ms time window before evidence onset
% It is centered on -fftlen = -150
MBrBin = -fftlen;
% Define variables that will be used in mixed linear effects model
MBatBL = []; RT = []; SAT = []; Contrast = []; Subject = []; IpsCon = []; LvsR = []; MBdiff = []; isCorr = [];

for subj = subjects
    clear allTrials thisMB allMB allRT
    allMB = nan(2,length(indicators.RT{subj})); allRT = nan(1,length(indicators.RT{subj}));
    
    % Only values that are not NaN can be entered when z-transforming motor
    % preparation amplitudes
    allTrials = find(~isnan(squeeze(nanmean(CSDstft{subj}(mueChans(lr),Frange,find(T==MBrBin),:),2))) & ~isnan(indicators.RT{subj})');
    
    % Within each subject and hemisphere, we z-transform all data points
    allMB(1,allTrials) = zscore(nanmean(CSDstft{subj}(mueChans(1),Frange,find(T==MBrBin),allTrials),2));
    allMB(2,allTrials) = zscore(nanmean(CSDstft{subj}(mueChans(2),Frange,find(T==MBrBin),allTrials),2));
    
    % Reaction times are z-scored within subjects
    allRT(1,allTrials) = zscore(indicators.RT{subj}(allTrials));
    
    for lr = 1:2
        
        trials = find(goodTrialsComb{subj}==1 &...
            squeeze(~isnan(nanmean(CSDstft{subj}(mueChans(lr),Frange,find(T==MBrBin),:),2)))' &...
            indicators.respLR{subj} == lr);
        
        for ic = 1:2
            
            if ic == 1
                thisMB = squeeze(allMB(3-lr,trials));
                thisDiff = squeeze(allMB(3-lr,trials)) - squeeze(allMB(lr,trials));
            else
                thisMB = squeeze(allMB(lr,trials));
                thisDiff = nan(1,length(trials));
            end
            RTS = allRT(1,trials);
            
            MBdiff = [MBdiff thisDiff];
            MBatBL = [MBatBL thisMB];
            RT = [RT RTS];
            SAT = [SAT satConds(indicators.cond{subj}(trials))];
            Subject = [Subject subj*ones(1,length(trials))];
            Contrast = [Contrast indicators.ContrLevels{subj}(trials)];
            
            isCorr = [isCorr indicators.respLR{subj}(trials) == indicators.LR{subj}(trials)];
            IpsCon = [IpsCon ic*ones(1,length(trials))];
            LvsR = [LvsR lr*ones(1,length(trials))];
        end
    end
end
save(fullfile(figData, 'MB_baselineLME'),'MBatBL','MBdiff','Subject','RT','SAT','Contrast','LvsR','isCorr','MBrBin','corrOnly')



%% Linear mixed-effects model for MB Excursion

% If not loaded already
% load(fullfile(figData, ['CSDstft_' conditions{1}]))
% load(fullfile(dataFolder, conditions{1}, ['CSDstftr_' conditions{1}]))

MBrBin = -fftlen;
MBpure = []; RT = []; SAT = []; Contrast = []; Subject = []; IpsCon = []; LvsR = [];

for subj = subjects
    
    clear allTrials allMB_Exc allRT
    allMB_Exc = nan(2,length(indicators.RT{subj}));
    allRT = nan(1,length(indicators.RT{subj}));
    allTrials = find((goodTrialsComb{subj}==1)' &...
        ~isnan(squeeze(nanmean(CSDstftr{subj}(mueChans(lr),Frange,find(Tr==MBrBin),:),2))) &...
        ~isnan(indicators.RT{subj})');
    allRT(1,allTrials) = zscore(indicators.RT{subj}(allTrials));
    
    for lr = 1:2
        thisConcat = zscore([squeeze(nanmean(CSDstftr{subj}(mueChans(lr),Frange,find(Tr==MBrBin),allTrials),2))'...
            squeeze(nanmean(CSDstft{subj}(mueChans(lr),Frange,find(T==MBrBin),allTrials),2))']);
        
        allMB_Exc(lr,allTrials) = thisConcat(1,1:length(allTrials))-thisConcat(1,length(allTrials)+1:end);
    end
    
    for lr = 1:2
        trials = find(goodTrialsComb{subj}==1 &...
            squeeze(~isnan(nanmean(CSDstftr{subj}(mueChans(lr),Frange,find(Tr==MBrBin),:),2)))' &...
            indicators.respLR{subj}==indicators.LR{subj} &...
            ~isnan(indicators.RT{subj}) &...
            indicators.respLR{subj} == lr);
        
        for ic = 1:2
            
            if ic == 1
                thisMB = squeeze(allMB_Exc(3-lr,trials));
            else
                thisMB = squeeze(allMB_Exc(lr,trials));
            end
            RTS = allRT(1,trials);
            
            MBpure = [MBpure thisMB];
            RT = [RT RTS];
            SAT = [SAT satConds(indicators.cond{subj}(trials))];
            Contrast = [Contrast indicators.ContrLevels{subj}(trials)];
            Subject = [Subject subj*ones(1,length(trials))];
            
            IpsCon = [IpsCon ic*ones(1,length(trials))];
            LvsR = [LvsR lr*ones(1,length(trials))];
        end
    end
end

ic = 2;
trials = find(IpsCon==ic);
datset=dataset({MBpure(trials)','MBexc'},...
    {Subject(trials)','Subject'},...
    {RT(trials)','RT'},...
    {SAT(trials)','SAT'},...
    {Contrast(trials)','Contrast'},...
    {LvsR(trials)','LvsR'});

predVar = 'MBexc';
fixedE = {'RT','RT:RT','SAT','LvsR','Contrast'};
randomE = {'Subject'};

save(fullfile(figData, 'ispiMuExc_LME'),'datset','predVar','fixedE','randomE');


%% Save data for JASP
trials = find(IpsCon==ic);
regression_table_MBexc = [[0 MBpure(trials)]' [1 RT(trials)]' [sqrt(2) RT(trials)]'.^2 [3 SAT(trials)]' [4 Contrast(trials)]'...
    [5 Subject(trials)]' [6 LvsR(trials)]'];
csvwrite(fullfile(figData, 'regression_table_MBexc.csv'), regression_table_MBexc)




















