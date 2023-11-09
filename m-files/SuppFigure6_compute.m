% SuppFigure6_compute

% Test for relationship between neural indicators and RT


% disp('Loading STFT SL...'); load(fullfile(figData, ['CSDstft_' conditions{1}]));
% disp('Loading STFT RL...'); load(fullfile(figData, ['CSDstftr_' conditions{1}]));
% disp('Loading CSDerpRIDE...'); load(fullfile(figData, ['CSDerpRIDE']))
% disp('Loading CSDerprRIDE...'); load(fullfile(figData, ['CSDerprRIDE']))


%% Mu/Beta and CPP and SSVEP: 3 RT-bins collapsed across conditions

nbins = 3;
load(fullfile(figData, ['CPP_MB_' num2str(nbins) 'RTbin']))
calc_newCPP = 0; calc_newCPPr = 0;
calc_newMB = 0; calc_newMBr = 0;
calc_newMBI = 1; calc_newMBrI = 1;
calc_SSVEP = 0; calc_SSVEPr = 0;


if calc_newCPP ==1;  bin2CPPslCs = nan(nbins,2,2,17,length(t)); end
if calc_newCPPr == 1; bin2CPPrlCs = nan(nbins,2,2,17,length(tr));end
if calc_newMB == 1; bin2MBslCs = nan(nbins,2,2,17,length(T));end
if calc_newMBr == 1; bin2MBrlCs = nan(nbins,2,2,17,length(Tr));end
if calc_newMBI == 1; bin2MBslCsI = nan(nbins,2,2,17,length(T));end
if calc_newMBrI == 1; bin2MBrlCsI = nan(nbins,2,2,17,length(Tr));end
if calc_SSVEP == 1; bin2SSVEPslCs = nan(nbins,2,2,17,length(T)); end
if calc_SSVEPr == 1; bin2SSVEPrlCs = nan(nbins,2,2,17,length(Tr)); end

meanRTs = nan(2,2,2,17);
for cc=1:2
    for l = 1:2
        
        thisERPrl = nan(max(subjects),3,2,length(tr),nbins);
        thisERPsl = nan(max(subjects),3,2,length(t),nbins);
        thisMBrl = nan(max(subjects),3,2,length(Tr),nbins);
        thisMBsl = nan(max(subjects),3,2,length(T),nbins);
        thisMBrlI = nan(max(subjects),3,2,length(Tr),nbins);
        thisMBslI = nan(max(subjects),3,2,length(T),nbins);
        thisSSVEPrl = nan(max(subjects),3,2,length(Tr),nbins);
        thisSSVEPsl = nan(max(subjects),3,2,length(T),nbins);

        thisRT = nan(max(subjects),3,2,nbins);
        for subj = subjects
            for d = 1:3
                for lr = 1:2
                    clear trlDelay rtlimDelay
                    
                    trlDelay = find(indicators.onsedelay{subj} == d &...
                        (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                        indicators.ContrLevels{subj} == l &...
                        indicators.LR{subj} == lr &...
                        indicators.LR{subj} == indicators.respLR{subj} &...
                        validrlockS{subj} &...
                        goodTrialsComb{subj});
                    
                    rtlimDelay = prctile(indicators.RT{subj}(trlDelay),[linspace(0,100,nbins+1)]);
                    for b=1:nbins
                        
                        condsCPPpooledBinned{cc,l,d,subj,lr,b} = find(indicators.RT{subj}>=rtlimDelay(b) &...
                            indicators.RT{subj}<rtlimDelay(b+1) &...
                            indicators.onsedelay{subj} == d &...
                            (indicators.cond{subj}==sats{cc}(1) | indicators.cond{subj}==sats{cc}(2)) &...
                            indicators.ContrLevels{subj} == l &...
                            indicators.LR{subj} == lr &...
                            indicators.LR{subj} == indicators.respLR{subj} &...
                            validrlockS{subj} &...
                            goodTrialsComb{subj});
                        
                        trials = [condsCPPpooledBinned{cc,l,d,subj,lr,b}];
                        
                        if calc_newCPPr == 1; thisERPrl(subj,d,lr,:,b) = squeeze(nanmean(CSDerprRIDE{subj}(chCPP,:,trials),3));end
                        if calc_newCPP ==1; thisERPsl(subj,d,lr,:,b) = squeeze(nanmean(CSDerpRIDE{subj}(chCPP,:,trials),3)); end
                        
                        if calc_newMBr == 1; thisMBrl(subj,d,lr,:,b) = squeeze(nanmean(nanmean(CSDstftr{subj}(mueChans(3-lr),Frange,:,trials),4),2));end
                        if calc_newMB == 1; thisMBsl(subj,d,lr,:,b) = squeeze(nanmean(nanmean(CSDstft{subj}(mueChans(3-lr),Frange,:,trials),4),2));end
                        
                        if calc_newMBrI == 1; thisMBrlI(subj,d,lr,:,b) = squeeze(nanmean(nanmean(CSDstftr{subj}(mueChans(lr),Frange,:,trials),4),2));end
                        if calc_newMBI == 1; thisMBslI(subj,d,lr,:,b) = squeeze(nanmean(nanmean(CSDstft{subj}(mueChans(lr),Frange,:,trials),4),2));end
                        
                        if calc_SSVEP == 1; thisSSVEPsl(subj,d,lr,:,b) = squeeze(nanmean(SSVEPsubtr{subj}(lr,:,trials)-SSVEPsubtr{subj}((3-lr),:,trials),3)); end
                        if calc_SSVEPr == 1; thisSSVEPrl(subj,d,lr,:,b) = squeeze(nanmean(SSVEPrSubtr{subj}(lr,:,trials)-SSVEPrSubtr{subj}((3-lr),:,trials),3)); end
                        
                        thisRT(subj,d,lr,b) = nanmean(indicators.RT{subj}(trials));
                    end
                end
            end
        end
        for b = 1:nbins
            
            if calc_newCPP ==1;  bin2CPPslCs(b,cc,l,:,:)  = squeeze(nanmean(nanmean(thisERPsl(:,:,:,:,b),2),3)); end
            if calc_newCPPr == 1; bin2CPPrlCs(b,cc,l,:,:)  = squeeze(nanmean(nanmean(thisERPrl(:,:,:,:,b),2),3));end
            
            if calc_newMB == 1; bin2MBslCs(b,cc,l,:,:)  = squeeze(nanmean(nanmean(thisMBsl(:,:,:,:,b),2),3)); end
            if calc_newMBr == 1; bin2MBrlCs(b,cc,l,:,:)  = squeeze(nanmean(nanmean(thisMBrl(:,:,:,:,b),2),3));end
            
            if calc_newMBI == 1; bin2MBslCsI(b,cc,l,:,:)  = squeeze(nanmean(nanmean(thisMBslI(:,:,:,:,b),2),3)); end
            if calc_newMBrI == 1; bin2MBrlCsI(b,cc,l,:,:)  = squeeze(nanmean(nanmean(thisMBrlI(:,:,:,:,b),2),3));end
            
            if calc_SSVEP == 1; bin2SSVEPslCs(b,cc,l,:,:)  = squeeze(nanmean(nanmean(thisSSVEPsl(:,:,:,:,b),2),3)); end
            if calc_SSVEPr == 1; bin2SSVEPrlCs(b,cc,l,:,:)  = squeeze(nanmean(nanmean(thisSSVEPrl(:,:,:,:,b),2),3)); end
            
            meanRTs(b,cc,l,:) = squeeze(nanmean(nanmean(thisRT(:,:,:,b),3),2));
        end
    end
end

save(fullfile(figData, ['CPP_MB+slope_' num2str(nbins) 'RTbin']),...
    'bin2CPPslCs','bin2CPPrlCs','bin2MBslCs','bin2MBrlCs','bin2SSVEPslCs','bin2SSVEPrlCs','meanRTs')

%% Compute significance for each time point

clear sigs p_RT
yy=1; f=1; p_RT{yy,f} = nan(1,length(T)); f=2; p_RT{yy,f} = nan(1,length(Tr));
yy=2; f=1; p_RT{yy,f} = nan(1,length(t)); f=2; p_RT{yy,f} = nan(1,length(tr));
yy=3; f=1; p_RT{yy,f} = nan(2,length(T)); f=2; p_RT{yy,f} = nan(2,length(Tr));


% SSVEP stimulus-locked
yy= 1; f=1;
sigs{yy,f} = nan(1,length(T));
for tpt = 1:length(T)
    
    psubj = nan(1,17);
    for subj = subjects
        amps = squeeze(nanmean(nanmean(nanmean(bin2SSVEPslCs(:,:,:,subj,tpt),2),3),5));
        rts = squeeze(nanmean(nanmean(meanRTs(:,:,:,subj),3),2));
        
        [R,P] = corrcoef(amps,rts);
        psubj(subj) = R(1,2);
        
    end
    [h,p,ci,stats] = ttest(psubj(subjects));
    sigs{yy,f}(1,tpt) = p;
    
end

% SSVEP response-locked
f=2;
sigs{yy,f} = nan(1,length(Tr));
for tpt = 1:length(Tr)
    
    psubj = nan(1,17);
    for subj = subjects
        amps = squeeze(nanmean(nanmean(nanmean(bin2SSVEPrlCs(:,:,:,subj,tpt),2),3),5));
        rts = squeeze(nanmean(nanmean(meanRTs(:,:,:,subj),3),2));
        
        [R,P] = corrcoef(amps,rts);
        psubj(subj) = R(1,2);
        
    end
    [h,p,ci,stats] = ttest(psubj(subjects));
    sigs{yy,f}(1,tpt) = p;
end

% CPP stimulus-locked
yy = 2; f=1;
sigs{yy,f} = nan(1,length(t)); binW = 50;
for tpt = 1:length(t)
    if tpt >= floor(binW/2) & tpt<length(t)-floor(binW/2)
        
        psubj = nan(1,17);
        for subj = subjects
            amps = squeeze(nanmean(nanmean(nanmean(bin2CPPslCs(:,:,:,subj,find(t>=t(tpt)-floor(binW/2) & t<=t(tpt)+floor(binW/2))),2),3),5));
            rts = squeeze(nanmean(nanmean(meanRTs(:,:,:,subj),3),2));
            
            [R,P] = corrcoef(amps,rts);
            psubj(subj) = R(1,2);
            
        end
        [h,p,ci,stats] = ttest(psubj(subjects));
        sigs{yy,f}(1,tpt) = p;
    end
end

% CPP response-locked
f=2;
sigs{yy,f} = nan(1,length(tr));
for tpt = 1:length(tr)
    if tpt >= floor(binW/2) & tpt<length(tr)-floor(binW/2)
        
        psubj = nan(1,17);
        for subj = subjects
            amps = squeeze(nanmean(nanmean(nanmean(bin2CPPrlCs(:,:,:,subj,find(tr>=tr(tpt)-floor(binW/2) & tr<=tr(tpt)+floor(binW/2))),2),3),5));
            rts = squeeze(nanmean(nanmean(meanRTs(:,:,:,subj),3),2));
            
            [R,P] = corrcoef(amps,rts);
            psubj(subj) = R(1,2);
            
        end
        [h,p,ci,stats] = ttest(psubj(subjects));
        sigs{yy,f}(1,tpt) = p;
        
    end
end


% MB stimulus-locked
yy=3; f=1;
sigs{yy,f} = nan(2,length(T));
for tpt = 1:length(T)
    
    for ic = 1:2
        
        psubj = nan(1,17);
        for subj = subjects
            if ic == 1
                amps = squeeze(nanmean(nanmean(nanmean(bin2MBslCs(:,:,:,subj,tpt),2),3),5));
            else
                amps = squeeze(nanmean(nanmean(nanmean(bin2MBslCsI(:,:,:,subj,tpt),2),3),5));
            end
            rts = squeeze(nanmean(nanmean(meanRTs(:,:,:,subj),3),2));
            
            [R,P] = corrcoef(amps,rts);
            psubj(subj) = R(1,2);
            
        end
        
        [h,p,ci,stats] = ttest(psubj(subjects));
        sigs{yy,f}(ic,tpt) = p;
        
    end
end

% MB response-locked
f=2;
sigs{yy,f} = nan(2,length(Tr));
for tpt = 1:length(Tr)
    
    for ic = 1:2
        
        psubj = nan(1,17);
        for subj = subjects
            if ic == 1
                amps = squeeze(nanmean(nanmean(nanmean(bin2MBrlCs(:,:,:,subj,tpt),2),3),5));
            else
                amps = squeeze(nanmean(nanmean(nanmean(bin2MBrlCsI(:,:,:,subj,tpt),2),3),5));
            end
            rts = squeeze(nanmean(nanmean(meanRTs(:,:,:,subj),3),2));
            
            [R,P] = corrcoef(amps,rts);
            psubj(subj) = R(1,2);
            
        end
        [h,p,ci,stats] = ttest(psubj(subjects));
        sigs{yy,f}(ic,tpt) = p;
    end
end

for yy = 1:3
    for f = 1:2
        p_RT{yy,f}(find(sigs{yy,f}<0.05)) = 1;
        
    end
end

save(fullfile(figData, [num2str(nbins) 'RTbin_signs']),'p_RT')

