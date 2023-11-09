% SuppFigure5_compute

% Mu/Beta amplitude just before evidence onset (window: -300 to 0ms
% before evidence onset) split into RT-bins plotted over RT


load(fullfile(figData, 'MB_BL_over_RT'), 'meanRTCtr', 'meanMu_BL', 'meanMu_Exc', 'meanMu_Exc')


% Plot traces with SEM separately for executed and withheld response
figure('Position', [100, 100, 1500, 400]); hold on
xlimBL = [400 1400];

clear dataMatrix dataMatrixExc dataMatrixRT dataMatrixRT_Exc
subjMean = squeeze(nanmean(nanmean(nanmean(nanmean(nanmean(nanmean(meanMu_BL,7),6),5),3),2),1));
subjMeanExc = squeeze(nanmean(nanmean(nanmean(nanmean(nanmean(meanMu_Exc,6),5),3),2),1));
subjMeanRT = squeeze(nanmean(nanmean(nanmean(nanmean(nanmean(meanRTCtr,6),5),3),2),1));
subjMeanRT_Exc = squeeze(nanmean(nanmean(nanmean(nanmean(nanmean(meanRT_Exc,6),5),3),2),1));
for subj = subjects
    dataMatrix(:,:,:,subj,:) = squeeze(nanmean(nanmean(meanMu_BL(:,:,:,subj,:,:,:),6),5)) - squeeze(repmat(subjMean(subj),CPPnumbins,2,2,2));
    dataMatrixExc(:,:,:,subj) = squeeze(nanmean(nanmean(meanMu_Exc(:,:,:,subj,:,:),6),5)) - squeeze(repmat(subjMeanExc(subj),CPPnumbins,2,2));
    
    dataMatrixRT(:,:,:,subj) = squeeze(nanmean(nanmean(meanRTCtr(:,:,:,subj,:,:),6),5)) - squeeze(repmat(subjMeanRT(subj),CPPnumbins,2,2));
    dataMatrixRT_Exc(:,:,:,subj) = squeeze(nanmean(nanmean(meanRT_Exc(:,:,:,subj,:,:),6),5)) - squeeze(repmat(subjMeanRT_Exc(subj),CPPnumbins,2,2));
end
for ic = 1:2
    subplot(1,3,ic); hold on
    for cc=1:2
        for l = 1:2
            
            plot(squeeze(nanmean(nanmean(nanmean(meanRTCtr(:,cc,l,:,:,:),6),5),4)),...
                squeeze(nanmean(nanmean(nanmean(meanMu_BL(:,cc,l,:,:,:,ic),6),5),4)),...
                'Color',colores{3-l,cc},'LineStyle',dash{ic},'Linewidth',3)
            
            for b=1:CPPnumbins
                
                v = nanstd(dataMatrix(b,cc,l,:,ic))./sqrt(length(subjects));
                m = nanmean(nanmean(nanmean(nanmean(meanMu_BL(b,cc,l,:,:,:,ic),6),5),4));
                
                vRT = nanstd(dataMatrixRT(b,cc,l,:))./sqrt(length(subjects));
                mRT = nanmean(nanmean(nanmean(nanmean(meanRTCtr(b,cc,l,:,:,:),6),5),4));
                
                line([mRT mRT],[m-v m+v],'Color',colores{3-l,cc},'LineStyle',dash{ic},'Linewidth',1)
                line([mRT-vRT mRT+vRT],[m m],'Color',colores{3-l,cc},'LineStyle',dash{ic},'Linewidth',1)
                
            end
        end
    end
    
    xlim(xlimBL); ylim([8 9.5]); axprefs(gca); set(gca,'Ydir','reverse');
    title(['Mu/Beta baseline plotted over RT']); xlabel('Time [ms]'); ylabel(['Mu/Beta (' choName{ic} ') [uV]'])
    
end

% Plot Mu/Beta Excursion for withheld response
subplot(1,3,3); hold on
for cc=1:2
    for l = 1:2
        
        plot(squeeze(nanmean(nanmean(nanmean(meanRT_Exc(:,cc,l,:,:,:),6),5),4)),...
            squeeze(nanmean(nanmean(nanmean(meanMu_Exc(:,cc,l,:,:,:),6),5),4)),...
            'Color',colores{3-l,cc},'LineStyle',dash{ic},'Linewidth',3)
        
        for b=1:CPPnumbins
            
            v = nanstd(dataMatrixExc(b,cc,l,:))./sqrt(length(subjects));
            m = nanmean(nanmean(nanmean(nanmean(meanMu_Exc(b,cc,l,:,:,:),6),5),4));
            
            vRT = nanstd(dataMatrixRT_Exc(b,cc,l))./sqrt(length(subjects));
            mRT = nanmean(nanmean(nanmean(nanmean(meanRT_Exc(b,cc,l,:,:,:),6),5),4));
            
            line([mRT mRT],[m-v m+v],'Color',colores{3-l,cc},'LineStyle',dash{ic},'Linewidth',1)
            line([mRT-vRT mRT+vRT],[m m],'Color',colores{3-l,cc},'LineStyle',dash{ic},'Linewidth',1)
            
        end
    end
end

xlim(xlimBL); ylim([-2 -.5]); axprefs(gca); set(gca,'Ydir','reverse');
title(['Mu/Beta excursion plotted over RT']); xlabel('Time [ms]'); ylabel(['Mu/Beta Excursion (' choName{ic} ') [uV]'])


