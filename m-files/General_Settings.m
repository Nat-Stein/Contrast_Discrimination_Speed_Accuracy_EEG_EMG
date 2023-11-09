% SATcontrast_GeneralSettings

% Add Matlab packages to path
addpath(fileFolder)
addpath(fullfile(matFolder, 'teg_repeated_measures_ANOVA')) 
addpath(fullfile(matFolder, 'sphertest')) 

% Number of subjects
nSubj = 17;

% EEG parameters
chans2interp = [69:72];
headchans = [1:68 73:96];       % EEG channels without reference
headchanels = [1:68 73:97];     % EEG channels with reference
emgchans = [69:72];             % EMG channels (EEG electrodes)
fs = 500;                       % Sampling frequency
nchan = 93;                     % Number of channels
chanlocs = readlocs(fullfile(fileFolder, 'JBhead96_sym.loc')); % Electrode locations (eeglab needs to be in path)


eplim = {[-1000 3200],[-500 4000]};     % Stimulus- and cue-locked epochs
trs = -500:300;         % Response-locked epoches in sample points (with respect to click)
tr = trs*1000/fs;       % Response-locked epochs in ms
t = [eplim{1}(1):2:eplim{1}(2)];        % Time for stimulus-locked epoch
tcue = [eplim{2}(1):2:eplim{2}(2)];     % Time for cue-locked epoch

% Frequency analysis
fftlen = 150;           % length of FFT
F = [0:25]*fs/fftlen;   % Frequencies analysed in STFT
T = [-200:50:1500];     % Time points of stimulus-locked STFT (center of time windows)
Tr = [-800:50:400];     % Time points of response-locked STFT (center of time windows)

beepThr = 15;           % Maximal difference between evidence onset and beep - use 10?

% ERPs will be extracted aligned to evidence onset (stimLock) and to the
% Speed/Accuracy cue (cueLocked)
conditions = {'stimLock','cueLocked'};


% Load behavioural measures
disp('Loading indicators...');
for c = 1:2
    thisFile = fullfile(figData, ['indicators_' conditions{c} '.mat']);
    if exist(thisFile,'file') == 2
        if c == 1
            load(thisFile);
            indicators_stimLock = indicators;
            load(fullfile(figData, ['indicatorsEYE_' conditions{c} '.mat']));
            load(fullfile(figData, 'goodTrials'))
        else
            load(thisFile);
        end
    end
end

indicators = indicators_stimLock;
% load(fullfile(figData, 'indicators_stimLock'))

load(fullfile(figData, 'elab'))     % Electrode locations

% Set general experimental measures
emgThresh = -250;
subjects = [1:15 17];
minRT = 0.2;                % Minimal RT for a trial to enter behavioural analysis
histbins = [0.2:0.08:2.4];  % RT bins for plotting accuracy as a function of RT

% Set RTs for all trials with no response to NaN
for subj = 1:17; indicators.RT{subj}(find(indicators.RT{subj}==0)) = NaN; end

colores = {[0 0 0.5625],[1 0.7 0],[.75 0 0],[34/255 100/255 20/255];...
    [0 0.6875 1.0000],[1 0.9 0],[1 0.0625 0],[60/255 190/255 113/255];...
    [0 0 0.9],[0.9 0.6 0],[.6 0.1 0],[160/255 240/255 120/255]};

dash = {'-','--','-.'};
sats = {[4 1];[2 3]};
satCond  = [1 2 2 1];
DlSlo = [1 2 1 2];


% ----------------------------------------------------------
% Parameters for plotting
CorrWrong = {'pooled','pooled';'correct','error'};
AccSpd = {'Accuracy','Speed'};
LoHi = {'Low Contrast','High Contrast'};
LeftRight = {'Left','Right'};
DlSlo_name = {'Deadline','Slope'};
choName = {'executed','withheld'};
xlimRL = [-800 200];
xlimSL = [-200 1500];
xlimRT = [400 1400];
widthBL = [0 500];

% ----------------------------------------------------------
% CPP
chCPP = 25;                 % Channel to record CPP
t_RT = [-130 -72];          % Time window for CPP at response
ylimitsCPP = [-5 30];       % Plotting limits
CPPnumbins=5; % For plots of signals over RT, trials will be split into CPPnumbins RT bins
% ----------------------------------------------------------

% ----------------------------------------------------------
% Mu/Beta amplitude
doLR = [1:2]; % [1:2] = average across left and right trials; 1 = only left trials; 2 = only right trials
Frange = find(F>=8 & F<30);     % Frequency range to of Mu/Beta
chL = strmatch('C3',elab);      % Channel for left-hemisphere Mu/Beta
chR = strmatch('C4',elab);      % Channel for right-hemisphere Mu/Beta
mueChans = [chL chR];           % Channels for both hemispheres' Mu/Beta
ylimitsMB = [6.5 10];           % Plotting limits
MBrBin = - fftlen;              % Mu/Beta amplitude at response will be measures in the 300ms time window just before the click
muBin = fftlen;                 % Mu/Beta amplitude aligned to stimulus will be measures in the 300ms time window just after evidence onset
% ----------------------------------------------------------

% ----------------------------------------------------------
% SSVEP
fftlenSSVEP = 200;
Fssvep = [0:25]*fs/fftlenSSVEP; % Frequencies for SSVEP derivation
chSSVEP = strmatch('Oz',elab,'exact');
ssvep_freqs = [20 40; 25 50];   % SSVEP frequencies and their first harmonics
ylimitsSSVEP = [-3 3];          % Plotting limits
% Determine data bin in which SSVEP frequencies and harmonics can be found
take_freqs{1} = []; take_freqs{2} = [];
for sfreq = 1:2
    [~, f1] = min(abs(ssvep_freqs(1,sfreq)-Fssvep));
    take_freqs{1} = [take_freqs{1} f1];
    [~, f1] = min(abs(ssvep_freqs(2,sfreq)-Fssvep));
    take_freqs{2} = [take_freqs{2} f1];
end
% ----------------------------------------------------------

% ----------------------------------------------------------
% Pupil
pupilLim = [-1 1.5];
pupNorm = 200;
% ----------------------------------------------------------

% ----------------------------------------------------------
% EMG
ylimitsEMG = {[0 65],[2.5 4.5]};
binWidthEMG = [10 10];
histOnsets = [0 600]; 
ylimsOnsets = {[0 0.3];[0 0.01]};
t_EMG = {-50 -175}; % Measurement time frame for EMG amplitude at response
fftlenEMG = 50; Femg = [0:25]*fs/fftlenEMG;
T2 = -200:25:1500; Tr2 = -800:25:400;
takeFemg = 2:find(Femg==150);
% for executed and withheld response
% ----------------------------------------------------------

% ----------------------------------------------------------
% Unilateral motor potentials
LRPChans = mueChans; 
ylimitsLRP = [-10 20];
% ----------------------------------------------------------


