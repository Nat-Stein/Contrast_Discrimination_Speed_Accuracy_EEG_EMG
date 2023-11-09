% SuppFigure5_stats

%% MuBeta at baseline (executed) - model comparison starting with no random slopes
load(fullfile(figData, 'MB_baselineLME'))

ic = 1; % Motor preparation towards the eventually executed response
trials = find(IpsCon==ic);

% Contruct data table for model
datset=dataset({MBatBL(trials)','MBatBL'},...
    {Subject(trials)','Subject'},...
    {RT(trials)','RT'},...
    {SAT(trials)','SAT'},...
    {Contrast(trials)','Contrast'},...
    {LvsR(trials)','LvsR'},...
    {isCorr(trials)','isCorr'});

predVar = 'MBatBL';
fixedE = {'RT','RT:RT','SAT','LvsR','Contrast'};
randomE = {'Subject'};

computeLME(datset, predVar, fixedE, randomE)

% RT slope: chi^2(2) = 18.9716, p = 7.5923e-05
% RT:RT slope: chi^2(2) = 14.2344, p = 0.00081103
% SAT slope: chi^2(2) = 1.1866, p = 0.55249
% LvsR slope: chi^2(2) = 7.276e-12, p = 1
% Contrast slope: chi^2(2) = 7.276e-12, p = 1
% Final model: MBatBL~RT+RT:RT+SAT+LvsR+Contrast+(1+RT+RT:RT|Subject)
% RT: chi^2(1) = 10.0865, p = 0.0014936
% RT:RT: chi^2(1) = 5.898, p = 0.015158
% SAT: chi^2(1) = 11.9121, p = 0.00055769
% LvsR: chi^2(1) = 0.00011578, p = 0.99141
% Contrast: chi^2(1) = 3.1177, p = 0.077444


%% MuBeta at baseline (withheld response)

ic = 2; % 2 = Motor preparation towards the withheld response
trials = find(IpsCon==ic);

% Contruct data table for model
datset=dataset({MBatBL(trials)','MBatBL'},...
    {Subject(trials)','Subject'},...
    {RT(trials)','RT'},...
    {SAT(trials)','SAT'},...
    {Contrast(trials)','Contrast'},...
    {LvsR(trials)','LvsR'},...
    {isCorr(trials)','isCorr'});

predVar = 'MBatBL';
fixedE = {'RT','RT:RT','SAT','LvsR','Contrast'};
randomE = {'Subject'};

computeLME(datset, predVar, fixedE, randomE)

% RT slope: chi^2(2) = 7.7146, p = 0.021125
% RT:RT slope: chi^2(2) = 7.5009, p = 0.023508
% SAT slope: chi^2(2) = 2.2863, p = 0.31881
% LvsR slope: chi^2(2) = 0, p = 1
% Contrast slope: chi^2(2) = -7.276e-12, p = 1
% Final model: MBatBL~RT+RT:RT+SAT+LvsR+Contrast+(1+RT+RT:RT|Subject)
% RT: chi^2(1) = 3.9073, p = 0.048076
% RT:RT: chi^2(1) = 1.3132, p = 0.25181
% SAT: chi^2(1) = 16.2497, p = 5.5519e-05
% LvsR: chi^2(1) = 0.11401, p = 0.73562
% Contrast: chi^2(1) = 1.4166, p = 0.23396




%% MB difference at baseline - model comparison starting without random slopes

trials = find(IpsCon==1);
datset=dataset({MBdiff(trials)','MBdiff'},...
    {Subject(trials)','Subject'},...
    {RT(trials)','RT'},...
    {SAT(trials)','SAT'},...
    {Contrast(trials)','Contrast'},...
    {LvsR(trials)','LvsR'},...
    {isCorr(trials)','isCorr'});

predVar = 'MBdiff';
fixedE = {'RT','RT:RT','SAT','LvsR','Contrast'};
randomE = {'Subject'};

computeLME(datset, predVar, fixedE, randomE)

% RT slope: chi^2(2) = 1.7486, p = 0.41715
% RT:RT slope: chi^2(2) = 0.32813, p = 0.84869
% SAT slope: chi^2(2) = -7.276e-12, p = 1
% LvsR slope: chi^2(2) = -7.276e-12, p = 1
% Contrast slope: chi^2(2) = 0.020605, p = 0.98975
% Final model: MBdiff~RT+RT:RT+SAT+LvsR+Contrast+(1|Subject)
% RT: chi^2(1) = 4.6494, p = 0.031064
% RT:RT: chi^2(1) = 3.2335, p = 0.072148
% SAT: chi^2(1) = 0.30226, p = 0.58247
% LvsR: chi^2(1) = 0.039284, p = 0.84289
% Contrast: chi^2(1) = 0.20371, p = 0.65174



%% Linear mixed-effects model for MB Excursion
load(fullfile(figData, 'ispiMuExc_LME'))

computeLME(datset, predVar, fixedE, randomE)

% RT slope: chi^2(2) = 11.4564, p = 0.0032529
% RT:RT slope: chi^2(2) = 0.037789, p = 0.98128
% SAT slope: chi^2(2) = 11.2929, p = 0.0035301
% LvsR slope: chi^2(2) = 71.9066, p = 2.2204e-16
% Contrast slope: chi^2(2) = 1.9638, p = 0.3746
% Final model: MBexc~RT+RT:RT+SAT+LvsR+Contrast+(1+RT+SAT+LvsR|Subject)
% RT: chi^2(1) = 2.7597, p = 0.096666
% RT:RT: chi^2(1) = 0.81677, p = 0.36613
% SAT: chi^2(1) = 2.3053, p = 0.12893
% LvsR: chi^2(1) = 8.1149, p = 0.0043903
% Contrast: chi^2(1) = 2.3566, p = 0.12475









