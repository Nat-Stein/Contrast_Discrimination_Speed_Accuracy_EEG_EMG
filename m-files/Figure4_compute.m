% Figure4_compute


% Dependin on the computer's memory, it may not be possible to load all
% four variables at the same time

disp('Loading STFT SL...'); load(fullfile(figData, ['CSDstft_' conditions{1} '.mat']));
disp('Loading STFT RL...'); load(fullfile(figData, ['CSDstftr_' conditions{1} '.mat']));
disp('Loading CSDerpRIDE...'); load(fullfile(figData, 'CSDerpRIDE'))
disp('Loading CSDerprRIDE...'); load(fullfile(figData, 'CSDerprRIDE'))

CPPnumbins=5; % For plots of signals over RT, trials will be split into CPPnumbins RT bins

%% Compute mean Stimulus-locked CPP
meanSL_CPP = nan(2,2,length(t));
% Compute traces for all four conditions
for cc = 1:2
    for l = 1:2
        thisERP = NaN*ones(3,2,max(subjects),length(t));
        for d = 1:3
            for lr = 1:2
                for subj = subjects
                    
                    trials = find(indicators.onsedelay{subj} == d  &...
                        (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                        indicators.ContrLevels{subj} == l &...
                        indicators.LR{subj} == lr &...
                        validrlockS{subj} &...
                        goodTrialsComb{subj});
                    
                    thisERP(d,lr,subj,:) = nanmean(CSDerpRIDE{subj}(chCPP,:,trials),3);
                end
            end
        end
        meanSL_CPP(cc,l,:) = squeeze(nanmean(nanmean(nanmean(thisERP,3),2),1));
    end
end

%% --------------------------------------------------------------------
%% Response-locked CPP
meanRL_CPP=NaN*ones(2,2,length(tr));
meanRL_CPP_stats=NaN*ones(2,2,max(subjects));
% Compute traces for all four conditions
for cc = 1:2
    for l = 1:2
        thisERP = NaN*ones(3,2,17,length(tr));
        for d = 1:3
            for lr = 1:2
                for subj = subjects
                    trials = find(indicators.onsedelay{subj} == d  &...
                        (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                        indicators.ContrLevels{subj} == l &...
                        indicators.LR{subj} == lr &...
                        validrlockS{subj} &...
                        goodTrialsComb{subj});
                    thisERP(d,lr,subj,:) = nanmean(CSDerprRIDE{subj}(chCPP,:,trials),3);
                end
            end
        end
        meanRL_CPP_stats(cc,l,:) = squeeze(nanmean(nanmean(nanmean(thisERP(:,:,:,find(tr>=t_RT(1) & tr<=t_RT(2))),1),2),4));
        meanRL_CPP(cc,l,:) = squeeze(nanmean(nanmean(nanmean(thisERP,3),2),1));
        
    end
end
%% --------------------------------------------------------------------
%% CPP at decision commitment (-130 to -70ms before click) plotted over RT
% Compute traces for all four conditions
meanRTCtr = NaN*ones(CPPnumbins,2,2,17,3,2);
meanCPPCtr = NaN*ones(CPPnumbins,2,2,17,3,2);
for cc=1:2
    for l = 1:2
        for subj = subjects
            for d = 1:3
                for lr = 1:2
                    clear trlDelay rtlimDelay
                    trlDelay = find(indicators.onsedelay{subj} == d &...
                        (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                        indicators.ContrLevels{subj} == l &...
                        indicators.LR{subj} == lr & validrlockS{subj} &...
                        goodTrialsComb{subj});
                    rtlimDelay = prctile(indicators.RT{subj}(trlDelay),[linspace(0,100,CPPnumbins+1)]);
                    for b=1:CPPnumbins
                        condsCPPpooledBinned{cc,l,d,subj,lr,b} = find(indicators.RT{subj}>=rtlimDelay(b) &...
                            indicators.RT{subj}<rtlimDelay(b+1) & indicators.onsedelay{subj} == d &...
                            (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                            indicators.ContrLevels{subj} == l & indicators.LR{subj} == lr &...
                            validrlockS{subj} &...
                            goodTrialsComb{subj});
                        trials = [condsCPPpooledBinned{cc,l,d,subj,lr,b}];
                        meanRTCtr(b,cc,l,subj,d,lr) = nanmean(indicators.RT{subj}(1,trials)*1000,2)';
                        meanCPPCtr(b,cc,l,subj,d,lr) = squeeze(nanmean(nanmean(CSDerprRIDE{subj}(chCPP,find(tr>=t_RT(1) & tr<=t_RT(2)),trials),3),2));
                    end
                end
            end
        end
    end
end

%% --------------------------------------------------------------------
%% Compute CPP topography at response
avgtopo = nan(17,97);
for subj = subjects
    trials = find(goodTrialsComb{subj});
    avgtopo(subj,:) = squeeze(nanmean(nanmean(CSDerprRIDE{subj}(:,find(tr>=t_RT(1) & tr<=t_RT(2)),trials),3),2));
end


save(fullfile(figData, 'plot_Fig3_CPP'), 'meanSL_CPP','avgtopo','meanRL_CPP','meanRL_CPP_stats','meanRTCtr','meanCPPCtr')

%% --------------------------------------------------------------------
%% Mu/Beta power
meanSL_mu = NaN*ones(2,2,2,length(T));
% Compute traces for all four conditions
for cc = 1:2
    for l = 1:2
        for ic = 1:2
            thisERP = NaN*ones(3,2,17,length(T));
            for d = 1:3
                for lr = 1:2
                    for subj = subjects
                        
                        trials = find(indicators.onsedelay{subj} == d &...
                            (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                            indicators.ContrLevels{subj} == l &...
                            indicators.respLR{subj} == lr &...
                            validrlockS{subj} &...
                            goodTrialsComb{subj});
                        
                        if ic == 1
                            thisERP(d,lr,subj,:) = nanmean(nanmean(CSDstft{subj}(mueChans(3-lr),Frange,:,trials),4),2);
                        else
                            thisERP(d,lr,subj,:) = nanmean(nanmean(CSDstft{subj}(mueChans(lr),Frange,:,trials),4),2);
                        end
                    end
                end
            end
            meanSL_mu(cc,l,ic,:) = squeeze(nanmean(nanmean(nanmean(thisERP,3),2),1));
            %meanBL_mu(cc,l,ic) = squeeze(nanmean(nanmean(nanmean(thisERP(:,:,:,find(T==MBrBin)),3),2),1));
        end
    end
end


%% --------------------------------------------------------------------
%% Topography of Mu/Beta amplitude:
% Difference of signal at response and baseline: motor preparation for the
% withheld response plotted on the left hemisphere and motor preparation
% for the executed response plotted on the right hemisphere
clear convElecs
for e = 1:size(chanlocs,2)-1
    allTheta(e) = chanlocs(e).theta;
    allRad(e) = chanlocs(e).radius;
end
for e = 1:size(chanlocs,2)-1
    if round(allTheta(e))~= 180 & round(allTheta(e))~= 0
        Matches = find(round(allRad*10) == round(allRad(e)*10) & allTheta >= -allTheta(e)-1 & allTheta <= -allTheta(e)+1);
        convElecs(e) = Matches;
    else
        convElecs(e) = e;
    end
end
diffLeft= nan(17,97);
diffRight = nan(17,97);
for subj = subjects
    
    trialsL = find(goodTrialsComb{subj} &...
        indicators.respLR{subj}==1);
    trialsR = find(goodTrialsComb{subj} &...
        indicators.respLR{subj}==2);
    
    % Topography of left-response trials: difference between amplitude at
    % response and at baseline
    diffLeft(subj,:) = squeeze(nanmean(nanmean(CSDstftr{subj}(:,Frange,find(Tr==MBrBin),trialsL),4),2)) -...
        squeeze(nanmean(nanmean(CSDstft{subj}(:,Frange,find(T==MBrBin),trialsL),4),2));
    
    % Topography of right-response trials: horizontal flip of difference
    % between amplitude at response and at baseline (left electrodes on the
    % riht and vice versa)
    diffRight(subj,:) = squeeze(nanmean(nanmean(CSDstftr{subj}(convElecs,Frange,find(Tr==MBrBin),trialsR),4),2)) -...
        squeeze(nanmean(nanmean(CSDstft{subj}(convElecs,Frange,find(T==MBrBin),trialsR),4),2));
end

%% --------------------------------------------------------------------
%% RL mu+beta
meanRL_mu = NaN*ones(2,2,2,length(Tr)); meanRL_MB_stats = NaN*ones(2,2,2,max(subjects));
% Compute traces for all four conditions and executed and withheld response
for cc = 1:2
    for l = 1:2
        for ic = 1:2
            
            thisERP = NaN*ones(3,2,17,length(Tr));
            for d = 1:3
                for lr = 1:2
                    for subj = subjects
                        
                        trials = find(indicators.onsedelay{subj} == d &...
                            (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                            indicators.ContrLevels{subj} == l &...
                            indicators.respLR{subj} == lr &...
                            validrlockS{subj} &...
                            goodTrialsComb{subj});
                        
                        if ic == 1
                            thisERP(d,lr,subj,:) = nanmean(nanmean(CSDstftr{subj}(mueChans(3-lr),Frange,:,trials),4),2);
                        else
                            thisERP(d,lr,subj,:) = nanmean(nanmean(CSDstftr{subj}(mueChans(lr),Frange,:,trials),4),2);
                        end
                    end
                end
            end
            meanRL_MB_stats(cc,l,ic,:) = squeeze(nanmean(nanmean(thisERP(:,:,:,find(Tr==MBrBin)),1),2));
            meanRL_mu(cc,l,ic,:) = squeeze(nanmean(nanmean(nanmean(thisERP,3),2),1));
            
        end
    end
end

%% Mu/Beta amplitude at the time of decision commitment (window: -300 to 0ms
% before button click) split into RT-bins plotted over RT
% Define position and dimenstions of subplot
% f=4; plotStart = x1st+plotWidth*sum(xWidths(1:f))/totalX + (f-1)*plotDist;
% plotX = plotWidth*xWidths(f+1)/totalX; plotY = yPos(yy);
% subplot('Position',[plotStart plotY plotX yHeight]); hold on
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
                            indicators.RT{subj}<rtlimDelay(b+1) & indicators.onsedelay{subj} == d &...
                            indicators.cond{subj}==c &...
                            indicators.ContrLevels{subj} == l &...
                            indicators.respLR{subj} == lr & validrlockS{subj} & goodTrialsComb{subj});
                    end
                end
            end
        end
    end
end

meanMuCtr = NaN*ones(CPPnumbins,2,2,17,3,2,2); meanRTCtr = NaN*ones(CPPnumbins,2,2,17,3,2);
for cc=1:2
    for l = 1:2
        for b=1:CPPnumbins
            for subj = subjects
                for d = 1:3
                    for lr = doLR
                        trials = [condsMuPooledBinnedM1{sats{cc}(1),l,d,subj,lr,b} condsMuPooledBinnedM1{sats{cc}(2),l,d,subj,lr,b}];
                        meanRTCtr(b,cc,l,subj,d,lr) = nanmean(indicators.RT{subj}(1,trials)*1000,2)';
                        meanMuCtr(b,cc,l,subj,d,lr,1) = squeeze(nanmean(nanmean(CSDstftr{subj}(mueChans(3-lr),Frange,find(Tr==MBrBin),trials),4),2));
                        meanMuCtr(b,cc,l,subj,d,lr,2) = squeeze(nanmean(nanmean(CSDstftr{subj}(mueChans(lr),Frange,find(Tr==MBrBin),trials),4),2));
                    end
                end
            end
        end
    end
end


save(fullfile(figData, 'plot_Fig3_MB'), 'meanSL_mu','diffRight','diffLeft','meanRL_mu','meanRL_MB_stats','meanRTCtr','meanMuCtr')




