% SuppFigure3_stats

%% 1. CPP amplitude at response (with Auditory Evoked Potential)

load(fullfile(figData, 'CPP_at_RT_stats_AEP'))

datset=dataset({CPPaep','CPPaep'},...
    {Subject','Subject'},...
    {RT','RT'},...
    {SAT','SAT'},...
    {Contrast','Contrast'},...
    {LvsR','LvsR'});

save(fullfile(figData, 'CPP_at_RT_stats_AEP'),'datset','fixedE','randomE','t_RT')

predVar = 'CPPaep';
fixedE = {'RT','RT:RT','SAT','LvsR','Contrast'};
randomE = {'Subject'};

computeLME(datset, predVar, fixedE, randomE)

% RT slope: chi^2(2) = 27.6265, p = 1.0023e-06
% RT:RT slope: chi^2(2) = 2.3467, p = 0.30933
% SAT slope: chi^2(2) = 4.2091, p = 0.1219
% LvsR slope: chi^2(2) = 6.7743, p = 0.033805
% Contrast slope: chi^2(2) = 13.2788, p = 0.0013078
% Final model: CPPride~RT+RT:RT+SAT+LvsR+Contrast+(1+RT+LvsR+Contrast|Subject)
% RT: chi^2(1) = 23.2173, p = 1.4469e-06
% RT:RT: chi^2(1) = 1.4483, p = 0.22881
% SAT: chi^2(1) = 4.2236, p = 0.039866
% LvsR: chi^2(1) = 1.1753, p = 0.27832
% Contrast: chi^2(1) = 2.3389, p = 0.12618


%% -----------------------------------------------------------------
%% 2. CPP slope stats using ANOVAs
% Compute rate of rise in CPP amplitude (with Auditory Evoked Potential)
% Build-up of CPP scales with Contrast and is increased under speed pressure

load(fullfile(figData, 'CPPSlope_AEP'),'CPPSlope_AEP')
load(fullfile(figData, 'CPPslopeSaep'),'CPPslopeSaep','tslope')
%% Run ANOVA on CPP slope per condition

clear cppSlope_300_50; con = 0;
for cc=1:2
    for l = 1:2
        s=0; con = con + 1;
        
        for subj = subjects
            s=s+1;
            
            cppSlope_300_50(s,con) = squeeze(nanmean(nanmean(CPPSlope_AEP(subj,cc,l,:,:),5),4));
            
        end
    end
end
output = teg_repeated_measures_ANOVA(cppSlope_300_50, [2,2], {'SAT','contrast'})
disp('CPP Slope (with AEP):')
for f = 1:length(output.labels)
    disp([output.labels{f} ': F(1,' num2str(output.R(f,3)) ')=' num2str(output.R(f,1)) '; p=' num2str(output.R(f,4))])
end
% SAT: F(1,15)=11.6405; p=0.0038629
% contrast: F(1,15)=21.0001; p=0.00035918
% SAT x contrast: F(1,15)=0.15919; p=0.69552



