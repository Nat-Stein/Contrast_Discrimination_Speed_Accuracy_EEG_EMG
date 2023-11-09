function [lateCAF_anova] = CAF_slope_stats(indicators, cutRTs,subjects)
% CAF_slope_stats splits trials with RT > cutRTs into 8 RT bins, fits
% slopes to the data points within experimental conditions, computes
% 3-Way ANOVA on the slopes and prints the results

clear condsBehCAlate
sats = {[4 1];[2 3]};
percentilesLate = [linspace(0,100,8+1)];
meanRTlate = NaN*ones(2,17,length(percentilesLate),2,2);
meanCorrlate = NaN*ones(2,17,length(percentilesLate),2,2);

for cc = 1:2
    for l = 1:2
        for lr = 1:2
            for subj = subjects
                
                condsBehCAlate{cc,l,subj,1,lr} = find((indicators.cond{subj} == sats{cc}(1) | indicators.cond{subj} == sats{cc}(2)) &...
                    indicators.respLR{subj}==indicators.LR{subj} &...
                    indicators.ContrLevels{subj} == l &...
                    indicators.LR{subj} == lr &...
                    indicators.RT{subj} >= cutRTs(subj,cc,l));
                
                condsBehCAlate{cc,l,subj,2,lr} = find((indicators.cond{subj} == sats{cc}(1) | indicators.cond{subj} == sats{cc}(2)) &...
                    indicators.respLR{subj}~=indicators.LR{subj} &...
                    indicators.ContrLevels{subj} == l &...
                    indicators.LR{subj} == lr &...
                    indicators.RT{subj} >= cutRTs(subj,cc,l));
                
                thisTrials = [condsBehCAlate{cc,l,subj,1,lr} condsBehCAlate{cc,l,subj,2,lr}];
                numCorr = length(condsBehCAlate{cc,l,subj,1,lr});
                thisRT = indicators.RT{subj}(thisTrials);
                pc = [prctile(thisRT,percentilesLate)]; clear ts
                for p = 1:length(pc)-1;
                    ts{p} = find(thisRT>=pc(p) & thisRT<=pc(p+1));
                    meanRTlate(lr,subj,p,cc,l) = nanmean(thisRT(ts{p}));
                    meanCorrlate(lr,subj,p,cc,l) = length(find(ts{p}<=numCorr))./length(ts{p});
                end
            end
        end
    end
end

% Fit slopes to all late CAFs

con = 0; clear lateCAF_anova
for cc = 1:2
    for l = 1:2
        for lr = 1:2
            s=0; con = con+1;
            for subj = subjects
                s=s+1;
                
                thisAcc = squeeze(meanCorrlate(lr,subj,1:8,cc,l));
                
                thisRT = squeeze(meanRTlate(lr,subj,1:8,cc,l));
                
                if sum(thisAcc) == 8
                    lateCAF_anova(s,con) = 0;
                else
                    p = polyfit(thisRT,thisAcc,1);
                    
                    lateCAF_anova(s,con) = p(1);
                end
            end
        end
    end
end

% Compute ANOVA and print output
output = teg_repeated_measures_ANOVA(lateCAF_anova, [2,2,2], {'Speed/Accuracy','Contrast','Left/Right Target'});
for comp = 1:length(output.labels)
    disp([output.labels{comp} ': F(1,15)=' num2str(output.R(comp,1)) ', p=' num2str(output.R(comp,4))])
end
disp('-------------------------------------------------')


end

