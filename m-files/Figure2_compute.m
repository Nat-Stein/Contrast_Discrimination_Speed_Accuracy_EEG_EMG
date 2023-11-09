%% Figure2_compute
% Compute all SSVEP and pupil measures needed to replicate Figure 2 and
% compute statistical analyses


%% Load stimulus-locked and response-locked ERPs

disp('Loading stimulus-locked ERPs...')
load(fullfile(figData, ['CSDerp_' conditions{1} '.mat']))
disp('Loading response-locked ERPs...')
load(fullfile(figData, ['CSDerpr_' conditions{1} '.mat']))
disp('Done loading...')

%% Calculate SSVEP - all SSVEP are calculated on non-CSDed ERPs
disp('Computing SSVEPs...')

for freq = 1:2
    
    freqBin{freq} = find(Fssvep == ssvep_freqs(freq,1) | Fssvep == ssvep_freqs(freq,2));
    neighbours{freq} = [find(Fssvep == ssvep_freqs(freq,1))-1 find(Fssvep == ssvep_freqs(freq,1))+1 find(Fssvep == ssvep_freqs(freq,2))-1,...
        find(Fssvep == ssvep_freqs(freq,2))+1];
end


clear SSVEPtopo
for subj = subjects
    disp(['Subject ' num2str(subj) '...'])
    
    for n = 1 : size(indicators.RT{subj},2)
        
        % Stimulus-locked SSVEP traces
        % Compute Short-Time Fourier Transform at all time points T
        clear this_stft
        for tt=1:length(T)
            temp = abs(fft(CSDerp{subj}(chSSVEP,find(t>=T(tt),1)-fftlenSSVEP/2+[1:fftlenSSVEP],n),[],2))./(fftlenSSVEP/2);
            this_stft(:,tt) = temp(:,1:length(Fssvep));
        end
        
        % Save entire frequency spectrum for later use in Supp Figure 2
        allFreq{subj}(:,:,n) = this_stft;
        
        % Calculate SSVEP as the difference between amplitude at stimulus
        % frequency and its first harmonic and the power at their
        % respective neighbouring frequencies
        % harmonic
        for freq = 1:2
            SSVEPsubtr{subj}(freq,:,n) = nanmean(this_stft(freqBin{freq},:)) - nanmean(this_stft(neighbours{freq},:));
        end
        
        % Response-locked SSVEP traces
        % Compute Short-Time Fourier Transform at all time points Tr
        clear this_stft
        for tt=1:length(Tr)
            temp = abs(fft(CSDerpr{subj}(chSSVEP,find(tr>=Tr(tt),1)-fftlenSSVEP/2+[1:fftlenSSVEP],n),[],2))./(fftlenSSVEP/2);
            this_stft(:,tt) = temp(:,1:length(Fssvep));
        end
        for freq = 1:2
            SSVEPrSubtr{subj}(freq,:,n) = nanmean(this_stft(freqBin{freq},:)) - nanmean(this_stft(neighbours{freq},:));
        end
        
        % SSVEP topography at baseline
        clear this_stft temp
        t_startBL = -800; t_endBL = 0;
        for e = 1:97
            temp = abs(fft(CSDerp{subj}(e,find(t>=t_startBL & t<t_endBL),n),[],2))./((t_endBL-t_startBL)/2);
            this_stft(:,e) = temp(:,1:length(Fssvep));
        end
        for freq = 1:2
            SSVEPtopoSubtrBL{subj}(freq,:,n) = nanmean(this_stft(freqBin{freq},:))-...
                nanmean(this_stft(neighbours{freq},:));
        end
        
    end
end


%% Save values for plotting and further computations
disp('Saving SSVEPs...')
save(fullfile(figData, 'SSVEPs'),'SSVEPsubtr','SSVEPrSubtr','SSVEPtopoSubtrBL',...
    't','tr','T','Tr','-v7.3') % 't_end','t_start',
% SSVEPtopoSubtrBL - topography in Figure 2
% SSVEPsubtr and SSVEPrSubtr - stimulus- and response-locked 
% traces

save(fullfile(figData, 'SSVEP_freqSpectrum'), 'T', 'allFreq')
% allFreq - analysis of individual SSVEP frequencies 

%% Calculating traces for plots and stats

load(fullfile(figData, 'pupilDiameters_SL_RL'))
% load(fullfile(figData, 'SSVEPs'))

doOverRT = 0;
plotRW = 1;

t_targ = round(eplim{1}(1)):2400;
t_resp = tr(1):tr(end);

for subj = subjects
    BlenS(subj) = max(indicators.trialNum{subj}); % Length of experimental block
    indicators.trialCount{subj} = (indicators.Block{subj}-1)*BlenS(subj)+indicators.trialNum{subj};
    
    BlenS(subj) = max(indicatorsEYE.trialNum{subj}); % Length of experimental block
    indicatorsEYE.trialCount{subj} = (indicatorsEYE.Block{subj}-1)*BlenS(subj)+indicatorsEYE.trialNum{subj};
end


plotRW = 0; rw = 1; 
meanSL_SSVEP = NaN*ones(2,2,2,length(T),plotRW+1);
meanBL_SSVEP = NaN*ones(2,2,max(subjects),plotRW+1);
meanRL_SSVEP = NaN*ones(2,2,2,length(Tr),plotRW+1);
meanRL_SSVEPS = NaN*ones(2,2,2,plotRW+1,max(subjects),length(Tr));
meanSL_SSVEPS = NaN*ones(2,2,2,plotRW+1,max(subjects),length(T)); % used for stats

meanSL_PUPIL = nan(2,2,2,length(t_targ));
meanBL_PUPIL = nan(2,2,max(subjects));
meanRL_PUPIL = nan(2,2,2,length(t_resp));
meanRL_PUPILS = nan(2,2,2,length(t_resp),max(subjects));
meanSL_PUPILS = nan(2,2,2,length(t_targ),max(subjects));

subtr_nbrsSSVEP = 1;
muBin = find(T==-fftlen);

for cc = 1:2
    for l = 1:2
        thisSSVEP = NaN*ones(3,2,17,length(T));
        thisSSVEP_RL = NaN*ones(3,2,17,length(Tr));
        thisSSVEP_BL = NaN*ones(3,2,17);
        
        thisEYE = NaN*ones(3,2,17,length(t_targ));
        thisEYE_RL = NaN*ones(3,2,17,length(t_resp));
        thisEYE_BL = NaN*ones(3,2,17);
        for lr = 1:2
            for d = 1:3
                for subj = subjects
                    
                    trialsSSVEP = find(indicators.onsedelay{subj} == d  &...
                        (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                        indicators.ContrLevels{subj} == l &...
                        indicators.LR{subj} == lr &...
                        validrlockS{subj} &...
                        goodTrialsComb{subj});
                    
                    % SSVEP
                    if subtr_nbrsSSVEP == 1
                        thisSSVEP_BL(d,lr,subj) = nanmean(SSVEPsubtr{subj}(1,find(T==-200),trialsSSVEP)-SSVEPsubtr{subj}(2,find(T==-200),trialsSSVEP),3);
                        thisSSVEP(d,lr,subj,:) = nanmean(SSVEPsubtr{subj}(1,:,trialsSSVEP)-SSVEPsubtr{subj}(2,:,trialsSSVEP),3);
                        thisSSVEP_RL(d,lr,subj,:) = nanmean(SSVEPrSubtr{subj}(1,:,trialsSSVEP)-SSVEPrSubtr{subj}(2,:,trialsSSVEP),3);
                    else
                        thisSSVEP_BL(d,lr,subj) = nanmean(SSVEP{subj}(1,find(T==-200),trialsSSVEP)-SSVEP{subj}(2,find(T==-200),trialsSSVEP),3);
                        thisSSVEP(d,lr,subj,:) = nanmean(SSVEP{subj}(1,:,trialsSSVEP)-SSVEP{subj}(2,:,trialsSSVEP),3);
                        thisSSVEP_RL(d,lr,subj,:) = nanmean(SSVEPr{subj}(1,:,trialsSSVEP)-SSVEPr{subj}(2,:,trialsSSVEP),3);
                    end
                    
                    trialsEYE = [];
                    for tt = 1:length(trialsSSVEP)
                        if length(find(indicatorsEYE.trialCount{subj} == indicators.trialCount{subj}(trialsSSVEP(tt)))) > 0
                            trialsEYE = [trialsEYE find(indicatorsEYE.trialCount{subj} == indicators.trialCount{subj}(trialsSSVEP(tt)))];
                        end
                    end
                    
                    thisEYE(d,lr,subj,:) = nanmean(pupilDiam{subj}(trialsEYE,:),1);
                    thisEYE_RL(d,lr,subj,:) = nanmean(pupilDiamR{subj}(trialsEYE,:),1);
                    thisEYE_BL(d,lr,subj) = nanmean(pupilDiam{subj}(trialsEYE,find(t_targ==0)),1);
                    
                end
            end
            
        end
        meanBL_SSVEP(cc,l,:,rw) = squeeze(nanmean(nanmean(thisSSVEP_BL,1),2))';
        meanBL_PUPIL(cc,l,:,rw) = squeeze(nanmean(nanmean(thisEYE_BL,1),2));
        
        meanSL_SSVEP(cc,l,:,:,rw) = squeeze(nanmean(nanmean(thisSSVEP,3),1));
        meanSL_PUPIL(cc,l,rw,:) = squeeze(nanmean(nanmean(nanmean(thisEYE,3),2),1));
        
        meanSL_SSVEPS(cc,l,:,rw,:,:) = squeeze(nanmean(thisSSVEP,1));
        meanSL_PUPILS(cc,l,rw,:,:) = squeeze(nanmean(nanmean(thisEYE,2),1))';
        
        meanRL_SSVEP(cc,l,:,:,rw) = squeeze(nanmean(nanmean(thisSSVEP_RL,3),1));
        meanRL_SSVEPS(cc,l,:,rw,:,:) = squeeze(nanmean(thisSSVEP_RL,1));
        meanRL_PUPIL(cc,l,rw,:) = squeeze(nanmean(nanmean(nanmean(thisEYE_RL,3),2),1));
        meanRL_PUPILS(cc,l,rw,:,:) = squeeze(nanmean(nanmean(thisEYE_RL,2),1))';
        
    end
end


save(fullfile(figData, 'stats_RL_SSVEPPupil'),'meanRL_SSVEPS','meanRL_PUPILS','Tr','t_resp',...
    'meanBL_SSVEP','meanBL_PUPIL','meanSL_PUPIL','meanSL_SSVEP','meanSL_PUPILS','meanSL_SSVEPS','T','t_targ','subtr_nbrsSSVEP')

meanRL_SSVEPS_CW = meanRL_SSVEPS;
save(fullfile(figData, 'meanRL_SSVEPS_CW'),'meanRL_SSVEPS_CW')

%% - Compare SSVEP traces over different RT bins and/or Correct/Wrong trials
% Calculating traces for plots and stats
% load(fullfile(figData, conditions{c}, 'SSVEPs'))

numbins = 1; percentiles = [linspace(0,100,numbins+1)];
meanSL_SSVEPS_bins = NaN*ones(2,2,2,2,max(subjects),length(T),numbins);
meanRT_SSVEPrw = NaN*ones(2,2,2,2,max(subjects),numbins);

for rw = 1:2
    for cc = 1:2
        for l = 1:2
            thisSSVEP = NaN*ones(3,2,17,length(T),numbins);
            meanRT_avg = NaN*ones(3,2,17,numbins);
            for lr = doLR
                for subj = subjects
                    for d = 1:3
                        
                        if rw == 1
                            trialsSSVEP = find(indicators.onsedelay{subj} == d  &...
                                (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                                indicators.ContrLevels{subj} == l &...
                                indicators.LR{subj} == lr &...
                                indicators.LR{subj} == indicators.respLR{subj} &... % correct response
                                validrlockS{subj} & goodTrialsComb{subj});
                        else
                            trialsSSVEP = find(indicators.onsedelay{subj} == d  &...
                                (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                                indicators.ContrLevels{subj} == l &...
                                indicators.LR{subj} == lr &...
                                indicators.LR{subj} ~= indicators.respLR{subj} &... % error
                                validrlockS{subj} & goodTrialsComb{subj});
                        end
                        
                        thisTrials = trialsSSVEP;
                        thisRT = indicators.RT{subj}(thisTrials);
                        pc = [prctile(thisRT,percentiles)]; clear ts
                        for p = 1:length(pc)-1;
                            ts{p} = find(thisRT>=pc(p) & thisRT<=pc(p+1));
                            meanRT_avg(d,lr,subj,p) = nanmean(thisRT(ts{p}));
                            
                            % SSVEP
                            thisSSVEP(d,lr,subj,:,p) = nanmean(SSVEPsubtr{subj}(1,:,thisTrials(ts{p}))-SSVEPsubtr{subj}(2,:,thisTrials(ts{p})),3);
                            
                        end
                        
                    end
                    for p=1:numbins
                        meanSL_SSVEPS_bins(cc,l,lr,rw,subj,:,p) = squeeze(nanmean(thisSSVEP(:,lr,subj,:,p),1));
                        meanRT_SSVEPrw(cc,l,lr,rw,subj,p) = squeeze(nanmean(meanRT_avg(:,lr,subj,p),1));
                    end
                end
                
            end
            
        end
    end
end

% Creates extra figure in Figure2_plot that is not shown in paper
save(fullfile(figData, 'SSVEP_RW'),'meanRT_SSVEPrw','meanSL_SSVEPS_bins','numbins')



