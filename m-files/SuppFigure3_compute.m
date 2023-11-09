% SuppFigure3_compute

% Statistics on the CPP amplitude at RT are computed in Figure4_stats

%% Load matrices needed for plotting all traces

% Dependin on the computer's memory, it may not be possible to load both
% variables at the same time

% disp('Loading CSDerp...'); load(fullfile(figData, ['CSDerp_' conditions{1}]))
% disp('Loading CSDerpr...'); load(fullfile(figData, ['CSDerpr_' conditions{1}]))

%% Plotting this

CPPnumbins=5; % For plots of signals over RT, trials will be split into 
% CPPnumbins RT bins

%%
meanSL_CPP = nan(2,2,length(t));
% Compute traces for all four conditions
for cc = 1:2
    for l = 1:2
        thisERP = NaN*ones(3,2,max(subjects),length(t));
        for d = 1:3
            disp([num2str(cc) ' - ' num2str(l) ' - ' num2str(d)])
            for lr = 1:2
                for subj = subjects
                    
                    trials = find(indicators.onsedelay{subj} == d  &...
                        (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                        indicators.ContrLevels{subj} == l &...
                        indicators.LR{subj} == lr &...
                        validrlockS{subj} &...
                        goodTrialsComb{subj});
                    
                    thisERP(d,lr,subj,:) = nanmean(CSDerp{subj}(chCPP,:,trials),3);
                    
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
disp('RL')
% Compute traces for all four conditions
for cc = 1:2
    for l = 1:2
        thisERP = NaN*ones(3,2,17,length(tr));
        for d = 1:3
            for lr = 1:2
                for subj = subjects
                    trials = [condsSSVEPpooled{sats{cc}(1),l,d,subj,lr} condsSSVEPpooled{sats{cc}(2),l,d,subj,lr}];
                    thisERP(d,lr,subj,:) = nanmean(CSDerpr{subj}(chCPP,:,trials),3);
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
disp('at RT')
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
                            indicators.ContrLevels{subj} == l &...
                            indicators.LR{subj} == lr &...
                            validrlockS{subj} &...
                            goodTrialsComb{subj});
                        trials = [condsCPPpooledBinned{cc,l,d,subj,lr,b}];
                        
                        meanRTCtr(b,cc,l,subj,d,lr) = nanmean(indicators.RT{subj}(1,trials)*1000,2)';
                        meanCPPCtr(b,cc,l,subj,d,lr) = squeeze(nanmean(nanmean(CSDerpr{subj}(chCPP,find(tr >= t_RT(1) & tr <= t_RT(2)), trials), 3), 2));
                    end
                end
            end
        end
    end
end

% --------------------------------------------------------------------
% Compute CPP topography at response
disp('topo')
avgtopo = nan(17,97);
for subj = subjects
    trials = find(goodTrialsComb{subj});
    avgtopo(subj,:) = squeeze(nanmean(nanmean(CSDerpr{subj}(:,find(tr >= t_RT(1) & tr <= t_RT(2)), trials), 3), 2));
end


save(fullfile(figData, 'plot_SuppFig3_CPPaep'), 'meanSL_CPP','avgtopo','meanRL_CPP','meanRL_CPP_stats','meanRTCtr','meanCPPCtr')





