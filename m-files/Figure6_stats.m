%% Figure6_stats



%% EMG slope

%% Measure slope of EMG muscle tone
load(fullfile(figData, ['emgR_' conditions{1}]))

% Compute EMG amplitude on a denser time scale
% STFT all EMG trials
Tr2_dense = -800:5:400;
clear emgSpecR_dense
for subj = subjects
    for n = 1:size(emgR{subj},3)
        for lr = 1:2            
            % Response-locked
            clear this_stft
            for tt=1:length(Tr2_dense)
                temp = abs(fft(emgR{subj}(lr,find(tr>=Tr2_dense(tt),1)-fftlenEMG/2+[1:fftlenEMG],n),[],2))./(fftlenEMG/2);
                this_stft(:,tt) = temp(:,1:length(Femg));
            end
            EMGr_dense{subj}(lr,n,:,:) = this_stft;
        end
    end
end

save(fullfile(figData, 'EMGr_dense'),'EMGr_dense','Tr2_dense')
%% Prepare data table for linear mixed-effects model

% load(fullfile(figData, conditions{1}, 'EMGr_dense'))

emgSlo = []; RT = []; SAT = []; Contrast = []; Subject = []; LvsR = [];
twin = [-175 -125]; tses = find(Tr2_dense>=twin(1) & Tr2_dense<=twin(2));


for subj = subjects
    allRTz{subj} = nan(size(indicators.RT{subj}));
    trials = ~isnan(indicators.RT{subj});
    allRTz{subj}(trials) = zscore(indicators.RT{subj}(trials));
end
for subj = subjects
    for cc = 1:2
        for l = 1:2
            
            for lr = 1:2
                
                trials = find((indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                    (~isnan(indicators.RT{subj})) &...
                    indicators.respLR{subj} == lr &...
                    goodTrialsComb{subj} &...
                    indicators.ContrLevels{subj} == l);
                
                for tt = 1:length(trials)
                    
                    % Measure Slope
                    thisHere = squeeze(nanmean(EMGr_dense{subj}(lr,trials(tt),takeFemg,tses),3));
                    pfit = polyfit(1:length(thisHere),thisHere',1);
                    
                    emgSlo = [emgSlo pfit(1)];
                    RT = [RT allRTz{subj}(trials(tt))];
                    
                    SAT = [SAT cc];
                    Contrast = [Contrast l];
                    Subject = [Subject subj];
                    LvsR = [LvsR lr];
                end
            end
        end
    end
end

datset=dataset({emgSlo','emgSlo'},...
    {Subject','Subject'},...
    {RT','RT'},...
    {SAT','SAT'},...
    {LvsR','LvsR'},...
    {Contrast','Contrast'});


predVar = 'emgSlo';
fixedE = {'RT','RT:RT','SAT','LvsR','Contrast'};
randomE = {'Subject'};

save(fullfile(figData, conditions{1}, 'EMG_slope_stats'),'datset','predVar','fixedE','randomE')

computeLME(datset, predVar, fixedE, randomE);

% RT slope: chi^2(2) = 109.3106, p = 0
% RT:RT slope: chi^2(2) = 37.9701, p = 5.6872e-09
% SAT slope: chi^2(2) = 92.2882, p = 0
% LvsR slope: chi^2(2) = 1170.5804, p = 0
% Contrast slope: chi^2(2) = 78.52, p = 0
% Final model: emgSlo~RT+RT:RT+SAT+LvsR+Contrast+(1+RT+RT:RT+SAT+LvsR+Contrast|Subject)
% RT: chi^2(1) = 10.0989, p = 0.0014836
% RT:RT: chi^2(1) = 8.3818, p = 0.0037899
% SAT: chi^2(1) = 8.236, p = 0.0041068
% LvsR: chi^2(1) = 0.090291, p = 0.76381
% Contrast: chi^2(1) = 0.86116, p = 0.35342


%% Linear mixed effects model of EMG amplitude at rersponse
load(fullfile(figData, ['EMGstftR_' conditions{1}]))

wEMG = []; RT = []; SAT = []; Contrast = []; Subject = []; Choice = []; LvsR = [];

satConds = [1 2 2 1];

emgBins = [-100 0;-175 -170];

for subj = subjects
    allRTz{subj} = nan(size(indicators.RT{subj}));
    trials = ~isnan(indicators.RT{subj});
    allRTz{subj}(trials) = zscore(indicators.RT{subj}(trials));
    for lr = 1:2
        for cho = 1:2
            
            clear trials
            if cho==1
                trials = find(indicators.respLR{subj} == lr & goodTrialsComb{subj});
            else
                trials = find(indicators.respLR{subj} == 3-lr & goodTrialsComb{subj});
            end
            
            thisEMG = squeeze(nanmean(nanmean(emgSpecR{subj}(lr,trials,takeFemg,find(Tr2>=emgBins(cho,1) & Tr2<emgBins(cho,2))),4),3));
            
            wEMG = [wEMG thisEMG];
            
            RT = [RT allRTz{subj}(trials)];
            SAT = [SAT satConds(indicators.cond{subj}(trials))];
            Contrast = [Contrast indicators.ContrLevels{subj}(trials)];
            Subject = [Subject subj*ones(1,length(trials))];
            Choice = [Choice cho*ones(1,length(trials))];
            LvsR = [LvsR lr*ones(1,length(trials))];
            
        end
    end
end


for cho = 1:2
    datset=dataset({wEMG(1,find(Choice==cho))','wEMG'},...
        {Subject(1,find(Choice==cho))','Subject'},...
        {RT(1,find(Choice==cho))','RT'},...
        {SAT(1,find(Choice==cho))','SAT'},...
        {LvsR(1,find(Choice==cho))','LvsR'},...
        {Contrast(1,find(Choice==cho))','Contrast'});
    
    predVar = 'wEMG';
    fixedE = {'RT','RT:RT','SAT','LvsR','Contrast'};
    randomE = {'Subject'};
    
    disp(['EMG at response - ' choName{cho}])
    
    computeLME(datset, predVar, fixedE, randomE)
    
end

% EMG at response - executed
% RT slope: chi^2(2) = 62.9727, p = 2.1205e-14
% RT:RT slope: chi^2(2) = 7.4481, p = 0.024136
% SAT slope: chi^2(2) = 70.6195, p = 4.4409e-16
% LvsR slope: chi^2(2) = 2775.0644, p = 0
% Contrast slope: chi^2(2) = 76.6101, p = 0
% Final model: wEMG~RT+RT:RT+SAT+LvsR+Contrast+(1+RT+RT:RT+SAT+LvsR+Contrast|Subject)
% RT: chi^2(1) = 5.418, p = 0.01993
% RT:RT: chi^2(1) = 11.3106, p = 0.00077064
% SAT: chi^2(1) = 18.1388, p = 2.0537e-05
% LvsR: chi^2(1) = 0.086598, p = 0.76855
% Contrast: chi^2(1) = 2.0368, p = 0.15354

% EMG at response - withheld
% RT slope: chi^2(2) = 18.2098, p = 0.00011112
% RT:RT slope: chi^2(2) = 14.9727, p = 0.00056069
% SAT slope: chi^2(2) = 7.8008, p = 0.020233
% LvsR slope: chi^2(2) = 245.5673, p = 0
% Contrast slope: chi^2(2) = 3.9603, p = 0.13805
% Final model: wEMG~RT+RT:RT+SAT+LvsR+Contrast+(1+RT+RT:RT+SAT+LvsR|Subject)
% RT: chi^2(1) = 8.6549, p = 0.0032619
% RT:RT: chi^2(1) = 1.8443, p = 0.17445
% SAT: chi^2(1) = 8.1066, p = 0.0044104
% LvsR: chi^2(1) = 0.22114, p = 0.63817
% Contrast: chi^2(1) = 0.41704, p = 0.51842


%% Mean delay between the last EMG onset burst before the button click and 
% the button click itself

load(fullfile(figData, 'EMG_plotMat'))

disp(['Mean MT (accuracy) = ' num2str(squeeze(nanmean(emgOn_subjMean(1,1,:),3))) ' +- ' num2str(squeeze(nanstd(emgOn_subjMean(1,1,:),[],3)))])
disp(['Mean MT (speed) = ' num2str(squeeze(nanmean(emgOn_subjMean(1,2,:),3))) ' +- ' num2str(squeeze(nanstd(emgOn_subjMean(1,2,:),[],3)))])
disp(['Mean MT (overall) = ' num2str(squeeze(nanmean(nanmean(emgOn_subjMean(1,:,:),3),2))) ' +- ' num2str(squeeze(nanstd(nanmean(emgOn_subjMean(1,:,:),2),[],3)))])
% Mean MT (accuracy) = -139.9247 +- 17.7471
% Mean MT (speed) = -130.2106 +- 16.7344
% Mean MT (overall) = -135.0676 +- 16.9506

choName = {'executed','withheld'};
cho = 1;
[h,p,ci,stats]  = ttest(emgOn_subjMean(cho,1,subjects),emgOn_subjMean(cho,2,subjects));
disp(['Onset times, ' choName{cho} ': t=' num2str(stats.tstat) ', p=' num2str(p)])
disp(['Accuracy - Speed = ' num2str(nanmean(emgOn_subjMean(cho,1,subjects)-emgOn_subjMean(cho,2,subjects)))...
    ' +- ' num2str(nanstd(emgOn_subjMean(cho,1,subjects)-emgOn_subjMean(cho,2,subjects)))])
% Onset times, executed: t=-6.0898, p=2.0696e-05
% Accuracy - Speed = -9.7141 +- 6.3806


%% Compute CPP peak time

% If not loaded already
% load(fullfile(figData, conditions{1}, 'CSDerprRIDE'))

smoothSpan = 101;
for subj = subjects
    for cc = 1:2
        for l = 1:2
            thisERP = NaN*ones(3,2,length(tr));
            for d = 1:3
                for lr = 1:2

                    trials = find(indicators.onsedelay{subj} == d  &...
                        (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                        indicators.ContrLevels{subj} == l &...
                        indicators.LR{subj} == lr &...
                        validrlockS{subj} &...
                        goodTrialsComb{subj});                
                    
                    thisERP(d,lr,:) = nanmean(CSDerprRIDE{subj}(chCPP,:,trials),3);
                end
            end
            meanSmooth = smooth(squeeze(nanmean(nanmean(thisERP,2),1)),smoothSpan,'lowess');
            allCPPmax(subj,cc,l) = max(meanSmooth(find(tr >= -150 & tr <= 100)));
            allTmax(subj,cc,l) = tr(find(meanSmooth == allCPPmax(subj,cc,l)));

        end
    end

end

save(fullfile(figData, 'CPP_peaks'),'allCPPmax','allTmax')

for cc = 1:2
    disp([AccSpd{cc} ': mean peak time +- sd = ' num2str(squeeze(nanmean(nanmean(allTmax(:,cc,:),3),1)))...
        ' +- ' num2str(squeeze(nanstd(nanmean(allTmax(:,cc,:),3),[],1)))])
end
% Accuracy: mean peak time +- sd = -44.2353 +- 38.4294
% Speed: mean peak time +- sd = -6.7647 +- 37.8839


%% 2-Way ANOVA shows that CPP peak is significantly closer to button click 
% under speed pressure.

clear cppPeak_anova; con = 0;
for cc=1:2
    for l = 1:2
        s=0; con = con+1;
        for subj = subjects
            s=s+1;
            cppPeak_anova(s,con) = squeeze(allTmax(subj,cc,l));
        end
    end
end
disp('CPP peak time')
output = teg_repeated_measures_ANOVA(cppPeak_anova, [2,2], {'SAT','contrast'})
for f = 1:length(output.labels)
   disp([output.labels{f} ': F(1,15)=' num2str(output.R(f,1)) '; p=' num2str(output.R(f,4))]) 
end

% SAT: F(1,15)=14.0391; p=0.0019433
% contrast: F(1,15)=0.31555; p=0.58259
% SAT x contrast: F(1,15)=1.8479; p=0.19411

%% Compute the difference between motor time and CPP peak time
load(fullfile(figData, 'emgOn_lastCorr'))

tDiff = nan(17,2,2);meanMT = nan(17,2,2); allTmaxMT = nan(17,2,2); 
smoothSpan = 101;
for subj = subjects
    for cc = 1:2
        for l = 1:2
            thisERP = NaN*ones(3,2,length(tr)); 
            thisMT = NaN*ones(3,2); 
            for d = 1:3
                for lr = 1:2
                    
                    % For this analysis, we are only using the trials which
                    % have a clearly visible EMG onset burst: ~isnan(lastCorr{subj})
                    trials = find(~isnan(lastCorr{subj}) &...
                        ~isnan(indicators.RT{subj}) &...
                        (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                        indicators.ContrLevels{subj} == l &...
                        indicators.onsedelay{subj} == d &...
                        indicators.respLR{subj} == lr &...
                        goodTrialsComb{subj});
                    
                    thisERP(d,lr,:) = nanmean(CSDerprRIDE{subj}(chCPP,:,trials),3);
                    thisMT(d,lr) = nanmean(lastCorr{subj}(trials));
                end
            end
            
            meanSmooth = smooth(squeeze(nanmean(nanmean(thisERP,2),1)),smoothSpan,'lowess');
            allCPPmaxMT(subj,cc,l) = max(meanSmooth(find(tr >=- 150 & tr <= 100))); 
            allTmaxMT(subj,cc,l) = tr(find(meanSmooth == allCPPmaxMT(subj,cc,l)));
            meanMT(subj,cc,l) = nanmean(nanmean(thisMT));
            
        end
    end
end

tDiff = allTmaxMT - meanMT;

clear cppPeakMT_anova; con = 0;
for cc=1:2
    for l = 1:2
        s=0; con = con+1;
        for subj = subjects
            s=s+1;
            cppPeakMT_anova(s,con) = squeeze(tDiff(subj,cc,l));
        end
    end
end
disp('Delay between EMG onset and CPP peak time:')
output = teg_repeated_measures_ANOVA(cppPeakMT_anova, [2,2], {'SAT','contrast'})
for f = 1:length(output.labels)
   disp([output.labels{f} ': F(1,15)=' num2str(output.R(f,1)) '; p=' num2str(output.R(f,4))]) 
end
% SAT: F(1,15)=7.8268; p=0.013524
% contrast: F(1,15)=0.0060636; p=0.93896
% SAT x contrast: F(1,15)=2.0651; p=0.17124


%% Number of EMG bursts in response withholding thumb

cho=2;
[h,p,ci,stats]  = ttest(emgOn_subjNum(cho,1,subjects),emgOn_subjNum(cho,2,subjects));
disp(['Num onsets, ' choName{cho} ': t = ' num2str(stats.tstat) ', p = ' num2str(p)])
disp(['# Accuracy - # Speed =' num2str(nanmean(emgOn_subjNum(cho,1,subjects)-emgOn_subjNum(cho,2,subjects)))...
    '+-' num2str(nanstd(emgOn_subjNum(cho,1,subjects)-emgOn_subjNum(cho,2,subjects)))])
disp(['# Speed = ' num2str(nanmean(emgOn_subjNum(cho,2,subjects))) ' +- ' num2str(nanstd(emgOn_subjNum(cho,2,subjects)))])
disp(['# Accuracy = ' num2str(nanmean(emgOn_subjNum(cho,1,subjects))) ' +- ' num2str(nanstd(emgOn_subjNum(cho,1,subjects)))])
% Num onsets, withheld: t = -4.6135, p = 0.00033785
% # Accuracy - # Speed =-14.125+-12.2468
% # Speed = 56.4375 +- 29.498
% # Accuracy = 42.3125 +- 23.6804


%% CPP at response is pocitive even for errors
load(fullfile(figData, 'CPPride_at_RT_stats'),'CPPride_noZ','Subject','RW')

errorCPP = nan(1,max(subjects));
for subj = subjects
    trials = find(Subject == subj & RW == 0);
    
    errorCPP(1,subj) = nanmean(CPPride_noZ(1,trials));
    
end
[h,p,ci,stats] = ttest(errorCPP);
disp(['t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p)])
% t(15)=3.63, p=0.0024
