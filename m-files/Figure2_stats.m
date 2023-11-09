% Figure2_stats

load(fullfile(figData, 'stats_RL_SSVEPPupil'),'meanSL_SSVEPS', 'meanRL_SSVEPS')

%% Stimulus-locked SSVEP 

disp('-------------------------')
% Stats stated in text for time window 250-450ms:
t_bin = [250 450];
clear ssvepSL_anova; con = 0;
for cc=1:2
    for rw = 1
        for l = 1:2
            for lr = 1:2
                s=0; con = con+1;
                for subj = subjects
                    s=s+1;
                    ssvepSL_anova(s,con) = squeeze(nanmean(meanSL_SSVEPS(cc,l,lr,rw,subj,find(T>= t_bin(1) & T<= t_bin(2))),6));
                end
            end
        end
    end
end
disp('SSVEP stim-locked: SAT, contrast, LvrR')
output = teg_repeated_measures_ANOVA(ssvepSL_anova, [2,2,2], {'SAT','contrast','LvsR'});

% Print results
for f = 1:size(output.R,1)
    disp([output.labels{f} ': F(1,15) = ' num2str(output.R(f,1)) ', p = ' num2str(output.R(f,4))])
end
% SAT: F(1,15) = 0.046732, p = 0.83176
% contrast: F(1,15) = 3.5236, p = 0.08009
% LvsR: F(1,15) = 26.7705, p = 0.0027646
% SAT x contrast: F(1,15) = 0.99077, p = 0.33534
% SAT x LvsR: F(1,15) = 5.4732, p = 0.033557
% contrast x LvsR: F(1,15) = 41.236, p = 1.15e-05
% SAT x contrast x LvsR: F(1,15) = 1.787, p = 0.20121


%% Test spericity of stimulus-locked SSVEP - this might not be working
allpsh = 0; clear allSpheres
[pWilk,pMauchly] = norm_sphere_test(ssvepSL_anova);
allpsh = allpsh+1; allSpheres.name{allpsh} = 'SSVEP Stimlocked 3-ANOVA'; allSpheres.pWilk{allpsh} = pWilk;
allSpheres.pMauchly{allpsh}=pMauchly; 
save(fullfile(figData, 'allSpheres_SSVEP'),'allSpheres')


%% Make big table of significance values for Stim-locked SSVEP power in different time windows
% Used for F-traces in Figure 3
load(fullfile(figData, 'stats_RL_SSVEPPupil'),'meanSL_SSVEPS')
clear allF
t_strt = -200:50:1500; winWidth = 0;
for tt = 1:length(t_strt)
    t_bin = [t_strt(tt) t_strt(tt)+winWidth];
    clear ssvepSL_anova; con = 0;
    for cc=1:2
        for rw = 1
            for l = 1:2
                for lr = 1:2
                    s=0; con = con+1;
                    for subj = subjects
                        s=s+1;
                        ssvepSL_anova(s,con) = squeeze(nanmean(meanSL_SSVEPS(cc,l,lr,rw,subj,find(T >= t_bin(1) & T <= t_bin(2))),6));
                    end
                end
            end
        end
    end
    output = teg_repeated_measures_ANOVA(ssvepSL_anova, [2,2,2], {'SAT','contrast','LvsR'});
    for i = 1:size(output.labels,2)
        allF{i}(tt,1) = output.R(i,1); 
        allF{i}(tt,2) = output.R(i,4);
    end
end

save(fullfile(figData, 'allF_SSVEP_SL'),'allF')

%% SSVEP baseline stats for the two frequencies separately

load(fullfile(figData, 'SSVEPs'))
load(fullfile(figData, 'goodTrials'))

disp('-------------------------')
SSVEP_BL12 = nan(17,2,2);
for freq = 1:2
    for subj = subjects
        for cc = 1:2
            trialsSSVEP = [];
            for lr = 1:2
                for d = 1:3
                    for sldl = 1:2
                        for l = 1:2
                            
                            addTrials = find(indicators.onsedelay{subj} == d  &...
                                indicators.cond{subj}==sats{cc}(sldl) &...
                                indicators.ContrLevels{subj} == l &...
                                indicators.LR{subj} == lr &...
                                validrlockS{subj} &...
                                goodTrialsComb{subj});

                            trialsSSVEP = [trialsSSVEP addTrials];
                        end
                    end
                end
            end
            SSVEP_BL12(subj,freq,cc) = squeeze(nanmean(SSVEPsubtr{subj}(freq,find(T==-200),trialsSSVEP)));
        end
    end
    [h,p,ci,stats] = ttest(squeeze(SSVEP_BL12(subjects,freq,1)),squeeze(SSVEP_BL12(subjects,freq,2)));
    disp(['SSVEP BL freq ' num2str(freq) ', Acc vs Spd: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p)])
    
end

% Results as in paper, varified April 1st, 2018:
% SSVEP BL freq 1, Acc vs Spd: t(15)=-0.42759, p=0.67502
% SSVEP BL freq 2, Acc vs Spd: t(15)=-0.52425, p=0.60776

% Save data for Bayes Factor analysis in JASP
ttest_table_bl20 = [[1:size(SSVEP_BL12,3)]; squeeze(SSVEP_BL12(subjects,1,:))];
csvwrite(fullfile(dataFolder, 'ttest_table_bl20.csv'),ttest_table_bl20)

ttest_table_bl25 = [[1:size(SSVEP_BL12,3)]; squeeze(SSVEP_BL12(subjects,2,:))];
csvwrite(fullfile(dataFolder, 'ttest_table_bl25.csv'),ttest_table_bl25)


%% Response-locked SSVEP 50ms before and after response time points around response
tpts = [-50 50];
disp('-------------------------')
for thisTr = 1:length(tpts)
    clear ssvepRL_anova; con = 0;
    for cc=1:2
        for rw = 1
            for l = 1:2
                for lr = 1:2
                    s=0; con = con+1;
                    for subj = subjects
                        s=s+1;
                        ssvepRL_anova(s,con) = squeeze(meanRL_SSVEPS(cc,l,lr,rw,subj,find(Tr==tpts(thisTr))));
                    end
                end
            end
        end
    end
    output = teg_repeated_measures_ANOVA(ssvepRL_anova, [2,2,2], {'SAT','contrast','LvsR'});
    % Print results
    disp('SSVEP resp-locked: SAT, contrast, LvrR')
    % thisTr = find(tpts==thisT);
    disp(['At ' num2str(tpts(thisTr))])
    for comp = 1:length(output.labels)
        disp([output.labels{comp} ': F=' num2str(output.R(comp,1)) ', p=' num2str(output.R(comp,4))]);
    end

end

% These are the stats that are in the paper
% At -50
% SAT: F=1.0666e-07, p=0.99974
% contrast: F=2.3162, p=0.14883
% LvsR: F=27.3834, p=0.002606
% SAT x contrast: F=0.014491, p=0.90578
% SAT x LvsR: F=5.6685, p=0.030961
% contrast x LvsR: F=38.0132, p=1.8075e-05
% SAT x contrast x LvsR: F=1.2794, p=0.27578
% At 50
% SAT: F=0.27907, p=0.60504
% contrast: F=1.7947, p=0.2003
% LvsR: F=26.4978, p=0.0082587
% SAT x contrast: F=4.6442, p=0.047816
% SAT x LvsR: F=2.8988, p=0.10927
% contrast x LvsR: F=41.4586, p=1.1158e-05
% SAT x contrast x LvsR: F=2.5781, p=0.1292

save(fullfile(figData, conditions{1}, 'SSVEP_stats_RL'), 'meanRL_SSVEPS') % , 'allF_resp'

%% Compute number of trials with RTs faster than 450ms

for subj = subjects
    rt_sub450(subj) = length(find(indicators.RT{subj} <= 450))./length(indicators.RT{subj});
end

disp('-------------------------')
disp(['mean+-std = ' num2str(mean(rt_sub450(subjects))) ' +- ' num2str(std(rt_sub450(subjects)))])
% mean+-std = 0.99866 +- 0.0025805
    
    
    
    
    