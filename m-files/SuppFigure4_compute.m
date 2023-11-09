% SuppFigure4_compute

% ------------------------------------------------------------------
%% Test robustness of effect of RT and Speed/Accuracy emphasis on CPP 
% amplitude at response for different time windows of measurement

% load(fullfile(figData, conditions{1}, 'CSDerprRIDE'))

t_stats = -600:10:200; wBin = 60;
clear p_FE chi2_FE
for tt = 1:length(t_stats)
    
    t_bin = [t_stats(tt) - wBin/2 t_stats(tt) + wBin/2];
    
    disp(['Computing: bin' num2str(tt) '/' num2str(length(t_stats))])
    
    CPPride = [];  RT = []; SAT = []; Contrast = []; Subject = []; LvsR = [];
    
    for subj = subjects
        
        trials = find(goodTrialsComb{subj}==1 &...
            squeeze(~isnan(nanmean(CSDerprRIDE{subj}(chCPP,find(tr >= t_bin(1) & tr < t_bin(2)),:),2)))');
        CPPsubj = zscore(squeeze(nanmean(CSDerprRIDE{subj}(chCPP,find(tr >= t_bin(1) & tr < t_bin(2)),trials),2)))';
        
        CPPride = [CPPride CPPsubj];
        
        RT = [RT zscore(indicators.RT{subj}(trials))];
        SAT = [SAT satConds(indicators.cond{subj}(trials))];
        Contrast = [Contrast indicators.ContrLevels{subj}(trials)];
        Subject = [Subject subj*ones(1,length(trials))];
        LvsR = [LvsR indicators.LR{subj}(trials)];
    end
    
    datset=dataset({CPPride','CPPride'},...
        {Subject','Subject'},...
        {RT','RT'},...
        {SAT','SAT'},...
        {Contrast','Contrast'},...
        {LvsR','LvsR'});
    
    finalMod = fitlme(datset,'CPPride~RT+RT:RT+SAT+LvsR+Contrast+(1+Contrast|Subject)');
    
    altMod = fitlme(datset,'CPPride~RT:RT+SAT+LvsR+Contrast+(1+Contrast|Subject)');
    results = compare(altMod, finalMod);
    p_FE(1,tt) = double(results(2,8));
    chi2_FE(1,tt) = double(results(2,6));
    
    altMod = fitlme(datset,'CPPride~RT+RT:RT+LvsR+Contrast+(1+Contrast|Subject)');
    results = compare(altMod, finalMod);
    p_FE(2,tt) = double(results(2,8));
    chi2_FE(2,tt) = double(results(2,6));
    
end

save(fullfile(figData, 'CPP_at_RT_twin'),'p_FE','chi2_FE','t_stats','wBin')

%% Compute response-locked unilateral motor potentials

meanRL_LRP = NaN*ones(2,2,2,length(tr));
for ic = 1:2
    for cc = 1:2
        for l = 1:2
            thisERP = NaN*ones(3,2,17,length(tr));
            for d = 1:3
                for lr = doLR
                    for subj = subjects
                        
                        trials = find(indicators.LR{subj} == indicators.respLR{subj} & indicators.onsedelay{subj} == d  &...
                            (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                            indicators.ContrLevels{subj} == l &...
                            indicators.respLR{subj} == lr & validrlockS{subj} & goodTrialsComb{subj});
                        
                        if ic == 1
                            thisERP(d,lr,subj,:) = nanmean(CSDerprRIDE{subj}(LRPChans(3-lr),:,trials),3);
                        else
                            thisERP(d,lr,subj,:) = nanmean(CSDerprRIDE{subj}(LRPChans(lr),:,trials),3);
                        end
                    end
                end
            end
            meanRL_LRP(cc,l,ic,:,rw) = squeeze(nanmean(nanmean(nanmean(thisERP,3),2),1));
        end
    end
end

save(fullfile(figData, 'meanRL_LRP'),'meanRL_LRP')

















































