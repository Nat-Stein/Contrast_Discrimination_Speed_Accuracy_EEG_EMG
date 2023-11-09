%% Figure4_stats



%% 1. MuBeta amplitude at response for executed response

load(fullfile(figData, 'MB_rtLME.mat'))

ic=1;
trials = find(IpsCon == ic);
datset=dataset({MBatRT(trials)','MBatRT'},...
    {Subject(trials)','Subject'},...
    {RT(trials)','RT'},...
    {SAT(trials)','SAT'},...
    {Contrast(trials)','Contrast'},...
    {LvsR(trials)','LvsR'});

predVar = 'MBatRT';
fixedE = {'RT','RT:RT','SAT','LvsR','Contrast'};
randomE = {'Subject'};

computeLME(datset, predVar, fixedE, randomE)
% RT slope: chi^2(2) = 35.6709, p = 1.7954e-08
% RT:RT slope: chi^2(2) = 22.5617, p = 1.2612e-05
% SAT slope: chi^2(2) = 18.1944, p = 0.00011198
% LvsR slope: chi^2(2) = 17.7878, p = 0.00013722
% Contrast slope: chi^2(2) = 25.1533, p = 3.4517e-06
% Final model: MBatRT~RT+RT:RT+SAT+LvsR+Contrast+(1+RT+RT:RT+SAT+LvsR+Contrast|Subject)
% RT: chi^2(1) = 0.42209, p = 0.5159
% RT:RT: chi^2(1) = 0.025524, p = 0.87307
% SAT: chi^2(1) = 0.54936, p = 0.45858
% LvsR: chi^2(1) = 0.71692, p = 0.39716
% Contrast: chi^2(1) = 1.9254, p = 0.16526



%% ------------------------------------------------------------------
%% 2. ANOVA for MuBeta amplitude at response: executed vs withheld response

load(fullfile(figData, 'MB_rtLME_executed_vs_withheld.mat'))

clear anovaInput; con = 0;
for cc=1:2
    for l = 1:2
        for ic = 1:2
            s=0; con = con+1;
            for subj = subjects
                s=s+1;
                anovaInput(s,con) = squeeze(meanMBresp(subj,cc,l,ic));
            end
        end
    end
end
output = teg_repeated_measures_ANOVA(anovaInput, [2,2,2], {'SAT','contrast','IpsCon'});
for f = 1:length(output.labels)
    disp([output.labels{f} ': F(1,15)=' num2str(output.R(f,1)) '; p=' num2str(output.R(f,4))])
end

% SAT: F(1,15)=3.078; p=0.099759
% contrast: F(1,15)=0.035778; p=0.85251
% IpsCon: F(1,15)=14.2593; p=0.0018304
% SAT x contrast: F(1,15)=0.25716; p=0.61945
% SAT x IpsCon: F(1,15)=0.00014296; p=0.99062
% contrast x IpsCon: F(1,15)=2.2975; p=0.15037
% SAT x contrast x IpsCon: F(1,15)=1.7491; p=0.2058

% ------------------------------------------------------------------
%% 3. CPP amplitude at baseline (before evidence onset)

load(fullfile(figData, 'CPP_blLME.mat'))

datset=dataset({CPPride','CPPride'},...
    {Subject','Subject'},...
    {RT','RT'},...
    {SAT','SAT'},...
    {Contrast','Contrast'},...
    {LvsR','LvsR'});

predVar = 'CPPride';
fixedE = {'RT','RT:RT','SAT','LvsR','Contrast'};
randomE = {'Subject'};

computeLME(datset, predVar, fixedE, randomE)

% RT slope: chi^2(2) = 4.3656e-11, p = 1
% RT:RT slope: chi^2(2) = 0.76724, p = 0.68139
% SAT slope: chi^2(2) = 1.6735e-10, p = 1
% LvsR slope: chi^2(2) = 0.91955, p = 0.63143
% Contrast slope: chi^2(2) = 2.3238, p = 0.31289
% Final model: CPPride~RT+RT:RT+SAT+LvsR+Contrast+(1|Subject)
% RT: chi^2(1) = 0.88991, p = 0.3455
% RT:RT: chi^2(1) = 0.011524, p = 0.91451
% SAT: chi^2(1) = 1.5148, p = 0.21841
% LvsR: chi^2(1) = 0.0013082, p = 0.97115
% Contrast: chi^2(1) = 0.041915, p = 0.83778


% ------------------------------------------------------------------
%% 4. CPP amplitude at response

calcWithAEP = 0;

if calcWithAEP == 0
    
    load(fullfile(figData, 'CPP_at_RT_stats'))
    
    predVar = 'CPPride';
    fixedE = {'RT','RT:RT','SAT','LvsR','Contrast'};
    randomE = {'Subject'};
    
    computeLME(datset, predVar, fixedE, randomE)
    
    % RT slope: chi^2(2) = 5.5095, p = 0.063626
    % RT:RT slope: chi^2(2) = 2.1275, p = 0.34515
    % SAT slope: chi^2(2) = 4.1289, p = 0.12689
    % LvsR slope: chi^2(2) = 5.0349, p = 0.080667
    % Contrast slope: chi^2(2) = 7.5236, p = 0.023242
    % Final model: CPPride~RT+RT:RT+SAT+LvsR+Contrast+(1+Contrast|Subject)
    % RT: chi^2(1) = 12.6748, p = 0.00037063
    % RT:RT: chi^2(1) = 1.4012, p = 0.23653
    % SAT: chi^2(1) = 3.1101, p = 0.077807
    % LvsR: chi^2(1) = 2.6024, p = 0.1067
    % Contrast: chi^2(1) = 2.3998, p = 0.12135
    
    % Save RIDE-corrected data for JASP analysis
    trials = 1:length(CPPride);
    regression_table_CPPrt = [[0 CPPride(trials)]' [1 RT(trials)]' [sqrt(2) RT(trials)]'.^2 [3 SAT(trials)]'...
        [4 Contrast(trials)]' [5 Subject(trials)]' [6 LvsR(trials)]'];
    csvwrite(fullfile(figData,'regression_table_CPPrt.csv'),regression_table_CPPrt)
    
else
    
    load(fullfile(figData, 'CPP_at_RT_stats_AEP'))
    
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
    
end



%% -------------------------------------------
% Compute rate of rise in CPP and Mu/Beta amplitude

%% 5. CPP slope stats using ANOVAs
% Build-up of CPP scales with Contrast and is increased under speed pressure

% Run ANOVA on CPP slope per condition

% If not loaded already
load(fullfile(figData, 'CPPSlopeRIDE'));

clear cppSlope_300_50; con = 0;
for cc=1:2
    for l = 1:2
        s=0; con = con + 1;
        
        for subj = subjects
            s=s+1;
            
            cppSlope_300_50(s,con) = squeeze(nanmean(nanmean(CPPSlopeRIDE(subj,cc,l,:,:),5),4));
            
        end
    end
end
disp('************************')
disp('CPP Slope: SAT, contrast')
output = teg_repeated_measures_ANOVA(cppSlope_300_50, [2,2], {'SAT','contrast'})
for f = 1:length(output.labels)
    disp([output.labels{f} ': F(1,' num2str(output.R(f,3)) ')=' num2str(output.R(f,1)) '; p=' num2str(output.R(f,4))])
end
disp('************************')
% SAT: F(1,15)=5.43; p=0.034165
% contrast: F(1,15)=12.9875; p=0.0026052
% SAT x contrast: F(1,15)=0.23353; p=0.6359



% ------------------------------------------------------------------
%% 6. Mu/Beta slope stats using ANOVAs

load(fullfile(figData, 'MBslopeS'));

%% Compute stats

ic = 1;
clear mbSlope_anova; con = 0;
for cc=1:2
    for l = 1:2
        s=0; con = con+1;
        for subj = subjects
            s=s+1;
            mbSlope_anova(s,con) = squeeze(nanmean(nanmean(MBSlope(subj,cc,l,:,:,ic),5),4));
        end
    end
end
disp('************************')
disp('MB Slope: SAT, contrast')
output = teg_repeated_measures_ANOVA(mbSlope_anova, [2,2], {'SAT','contrast'})
for f = 1:length(output.labels)
    disp([output.labels{f} ': F(' num2str(output.R(f,2)) ',' num2str(output.R(f,3)) ')=' num2str(output.R(f,1)) '; p=' num2str(output.R(f,4))])
end
disp('************************')
% SAT: F(1,15)=11.5067; p=0.0040221
% contrast: F(1,15)=9.2173; p=0.0083381
% SAT x contrast: F(1,15)=0.57975; p=0.45821










