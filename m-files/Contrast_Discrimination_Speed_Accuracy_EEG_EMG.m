% Contrast_Discrimination_Speed_Accuracy_EEG_EMG
% Copyright (c) 2023 Natalie Steinemann

% This file is a compilation of all files requires to computes all stats and plots all figures


%% Specify file locations
% Location of folder 'Contrast_Discrimination_Speed_Accuracy_EEG_EMG'
prefix_1 = '';
% Location of Matlab packages
matFolder = '\MATLAB\';
% Download toolbox norm_sphere_test for spericity testing: https://www.mathworks.com/matlabcentral/fileexchange/3694-sphertest

% This is the directory of the folder for both data and m-files
prefix_ = [prefix_1 'Contrast_Discrimination_Speed_Accuracy_EEG_EMG'];
fileFolder = fullfile(prefix_, 'm-files');
cd(fileFolder)

dataFolder = fullfile(prefix_, 'Data'); 
figData = fullfile(dataFolder, 'Data for Figures'); 

%% General settings across analyses
% General_Settings.m needs to be run before starting any data analysis

General_Settings

%% ---------------------------------------------------------
%% Statistics and plots by figure
%% ---------------------------------------------------------

%% Figure 1: Behaviour - histogrammes and Conditional Accuracy

Figure1_plot

% Compute statistics on behavioural measures
Behaviour_stats

%% Figure 2: SSVEP

% Figure2_compute

Figure2_stats

Figure2_plot

SSVEP_CorrVsErr % ANOVA testing behavioural effect of SSVEP


%% Figure 3: Pupil size

Figure3_stats
% 1. Pupil size at baseline
% 2. Effect of speed pressure on pupil size increases over time
% 3. Influence of pupil size on SSVEP

Figure3_plot
% Data matrix computed in Figure3_stats

%% Figure 4: CPP and Mu+Beta at baseline, stim-locked, response-locked and over RT
% Plot cumulative reaction time histogramme and conditional accuracy
% for speed and accuracy for high and low contrast

% Figure4_compute

% Figure4_prepstats

Figure4_stats
% 1. Mu/Beta amplitude at RT for executed response (LME)
% 2. Contralateral vs ipsilateral Mu/Beta amplitude at response
% 3. CPP amplitude at baseline (LME)
% 4. CPP amplitude at response (LME)
% - with and without auditory-evoked
% potential
% 5. CPP rate of rise (ANOVA)
% 6. Mu/Beta rate of rise (ANOVA)

Figure4_plot

CAF_tail_stats


%% Figure 6: EMG

% Figure6_compute 

Figure6_plot

Figure6_stats
% Includes analysis of delay between CPP peak time and button click



%% Supplementary Figure 1 - behaviour split be Deadline and Slope 
% conditions within each Speed/Accuracy Regime

SuppFigure1_plot

% Also includes Komolgorov-Smirnov to test whether RT distributions within 
% behavioural conditions are significantly different between Deadline
% and Slope conditions


%% Supplementary Figure 2

SuppFigure2_compute
% Also incluses stats

SuppFigure2_plot

%% Supplementary Figure 3

% SuppFigure3_compute 

SuppFigure3_plot

% SuppFigure3_prepstats

SuppFigure3_stats
% 1. CPP with AEP at RT (LME)
% 2. CPP with AEP rate of rise (ANOVA)


%% Supplementary Figure 4

% SuppFigure4_compute

SuppFigure4


%% Supplementary Figure 5

% SuppFigure5_compute

SuppFigure5_plot

% SuppFigure5_prepstats

SuppFigure5_stats

%% Supplementary Figure 6

% SuppFigure6_compute

SuppFigure6_plot



