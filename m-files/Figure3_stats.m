% Figure3_stats


%% 1. Pupil size at baseline
load(fullfile(figData, 'stats_RL_SsvepPupil'),'meanBL_PUPIL')

[h,p,ci,stats] = ttest(squeeze(nanmean(meanBL_PUPIL(1,:,:),2)),squeeze(nanmean(meanBL_PUPIL(2,:,:),2)));
disp(['Basline Pupil (correct) Acc vs Spd: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p)])
% Basline Pupil (correct) Acc vs Spd: t(15)=-2.245, p=0.040273


%% ------------------------------------------------------------
%% 2. Effect of speed pressure on pupil size increases over time
% Interaction between time and SAT

clear pupilT_anova
con=0;
for cc=1:2
    for tt = find(T>0)
        s=0; con = con+1;
        for subj = subjects
            s=s+1;
            pupilT_anova(s,con) = squeeze(nanmean(nanmean(meanSL_PUPILS(cc,:,1,find(t_targ>=T(tt)-25 & t_targ<T(tt)+25),subj),2),4));
        end
    end
end
disp('PUPIL: SAT, time')
output = teg_repeated_measures_ANOVA(pupilT_anova, [2,size(pupilT_anova,2)/2], {'SAT','Time'});
for f = 1:size(output.R,1)
    disp([output.labels{f} ': F(29,435) = ' num2str(output.R(f,1)) ', p = ' num2str(output.R(f,4))])
end
% [pWilk,pMauchly] = norm_sphere_test(pupilT_anova);
% allpsh = allpsh+1; allSpheres.name{allpsh}=['Pupil over time ANOVA']; allSpheres.pWilk{allpsh}=pWilk;
% allSpheres.pMauchly{allpsh}=pMauchly; save(fullfile(resultsFolder, 'allSpheres_SSVEP'),'allSpheres')

% SAT: F(29,435) = 23.5483, p = 0.030639
% Time: F(29,435) = 41.9185, p = 2.1724e-06
% SAT x Time: F(29,435) = 17.6185, p = 0.00045951


%% ----------------------------------------------------------
%% 3. Influence of pupil size on SSVEP
% Split trials into two bins median-split by pupil size
load(fullfile(figData, 'SSVEPs'),'SSVEPsubtr')
load(fullfile(figData, 'pupilDiameters_SL_RL'), 'pupilDiam')

numBins = 2;

t_targ = round(eplim{1}(1)):round(eplim{1}(2));

meanSL_SSVEP = NaN*ones(2,2,2,length(T));
meanSL_SSVEPbin = NaN*ones(2,2,2,2,length(T));
meanSL_SSVEPbinS = NaN*ones(2,2,2,2,17,length(T));

for cc = 1:2
    for l = 1:2
        thisSSVEP = NaN*ones(3,2,17,length(T));
        thisSSVEPbin = NaN*ones(2,3,2,17,length(T));
        for lr = doLR
            for d = 1:3
                for subj = subjects
                    % trialsSSVEP = [condsSSVEPpooled{sats{cc}(1),l,d,subj,lr} condsSSVEPpooled{sats{cc}(2),l,d,subj,lr}];
                    trialsSSVEP = [];
                    for ccc = 1 : 2 % Combine both methods of Speed/Accuracy encouragement
                        addTrials = find(indicators.onsedelay{subj} == d  &...
                        indicators.cond{subj}==sats{cc}(ccc) &...
                        indicators.ContrLevels{subj} == l &...
                        indicators.LR{subj} == lr &...
                        validrlockS{subj} &...
                        goodTrialsComb{subj});
                        trialsSSVEP = [trialsSSVEP addTrials];
                    end
                    
                    % Pupil trials
                    trialsEYE = [];
                    for tt = 1:length(trialsSSVEP)
                        if length(find(indicatorsEYE.trialCount{subj} == indicators.trialCount{subj}(trialsSSVEP(tt)))) > 0
                            trialsEYE = [trialsEYE find(indicatorsEYE.trialCount{subj} == indicators.trialCount{subj}(trialsSSVEP(tt)))];
                        end
                    end
                    
                    pupilArea = squeeze(nanmean(pupilDiam{subj}(trialsEYE,find(t_targ>=0 & t_targ<=1500)),2));
                    [pupD,ind] = sort(pupilArea);
                    ntr = length(pupilArea);
                    for b = 1:numBins
                        bTrials = ind([(b-1)*floor(ntr/2)+1:b*floor(ntr/2)]);
                        thisSSVEPbin(b,d,lr,subj,:) = nanmean(SSVEPsubtr{subj}(1,:,trialsSSVEP(bTrials))-SSVEPsubtr{subj}(2,:,trialsSSVEP(bTrials)),3);
                        
                    end
                    
                    % SSVEP
                    thisSSVEP(d,lr,subj,:) = nanmean(SSVEPsubtr{subj}(1,:,trialsSSVEP)-SSVEPsubtr{subj}(2,:,trialsSSVEP),3);
                    
                end
            end
        end
        meanSL_SSVEPbinS(:,cc,l,:,:,:) = squeeze(nanmean(thisSSVEPbin,2));
    end
end

save(fullfile(figData, 'SSVEPs_by_Pupil'),'meanSL_SSVEPbinS')

%% 3. Influence of pupil size on SSVEP
% Compute ANOVA on effect of pupil size on SSVEP

% load(fullfile(figData, 'SSVEPs_by_Pupil'),'meanSL_SSVEPbinS')

t_bin = [350 550];
clear ssvepSL_anova; con = 0;
for cc=1:2
    for l = 1:2
        for lr = 1:2
            for b = 1:2
                s=0; con = con+1;
                for subj = subjects
                    s=s+1;
                    ssvepSL_anova(s,con) = squeeze(nanmean(meanSL_SSVEPbinS(b,cc,l,lr,subj,find(T>= t_bin(1) & T<= t_bin(2))),6));
                end
            end
        end
    end
end
output = teg_repeated_measures_ANOVA(ssvepSL_anova, [2,2,2,2], {'SAT','contrast','LvsR','Pupil'});
disp(' ')
disp('*************************')
disp(['Influence of pupil size on SSVEP w. SAT, contrast, LvR - ' num2str(t_bin(1)) '-' num2str(t_bin(2)) 'ms'])
f = 10; disp([output.labels{f} ': F(1,15) = ' num2str(output.R(f,1)) ', p = ' num2str(output.R(f,4))])
disp('*************************')
for f = 1:size(output.R,1)
    disp([output.labels{f} ': F(1,15) = ' num2str(output.R(f,1)) ', p = ' num2str(output.R(f,4))])
end

% *************************
% Influence of pupil size on SSVEP w. SAT, contrast, LvR - 350-550ms
% LvsR x Pupil: F(1,15) = 10.8287, p = 0.0049525
% *************************
% SAT: F(1,15) = 0.037903, p = 0.84825
% contrast: F(1,15) = 2.7806, p = 0.11615
% LvsR: F(1,15) = 27.5657, p = 0.0025611
% Pupil: F(1,15) = 1.1204, p = 0.30657
% SAT x contrast: F(1,15) = 0.30792, p = 0.58714
% SAT x LvsR: F(1,15) = 10.7135, p = 0.0051339
% SAT x Pupil: F(1,15) = 0.091768, p = 0.7661
% contrast x LvsR: F(1,15) = 39.6761, p = 1.4265e-05
% contrast x Pupil: F(1,15) = 1.1745, p = 0.29558
% LvsR x Pupil: F(1,15) = 10.8287, p = 0.0049525
% SAT x contrast x LvsR: F(1,15) = 1.3173, p = 0.26905
% SAT x contrast x Pupil: F(1,15) = 0.12426, p = 0.72936
% SAT x LvsR x Pupil: F(1,15) = 3.1047, p = 0.098429
% contrast x LvsR x Pupil: F(1,15) = 1.2862, p = 0.27455
% SAT x contrast x LvsR x Pupil: F(1,15) = 1.7302, p = 0.20814

