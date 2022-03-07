clear all
close all
clc

%This is the code used to analyse the expression data and to fit it to the
%gene expression model

%Input of the code is a table (in this example 'meansCUP') that was
%extracted using YeaZ

%The structure of the table is as follow:
%A row contains expression values for a single cell through time
%row id is the id of a cell as defined in the YeaZ output
%column id is the timepoint (frame number) in the analyzed movie 
%(0 timepoint corresponds to colum index 1 since counting in Matlab starts from 1)

%For expression values before the birth of a cell,
%and in case the cell with a specific id doesn't exist
%the table contains NaN values

load meansCUP
t_start_means = 33+1; %T = -60 min corresponds to the frame 33; +1 because Matlab starts counting from 1
m_cup = means(:, t_start_means:end); %
clear means
means = m_cup;


frequency_of_imaging = 10; %frequency of imaging; in minutes
starting_frame = -6;
ending_frame = 39;
%defining time axis
time = starting_frame*frequency_of_imaging : frequency_of_imaging : ending_frame*frequency_of_imaging;
average = calculate_average(means);

frames_on = 2:10; %frames on which tau-on and i are extracted
frames_base = 1:7; %frames during which there is no activation (-60 to 0). used to calculate base
frames_off = 28:34; %frames for finding tau-off
frames_tail = 34:46; %frames for fitting degradation rate

plotting_average = 1;

x0 = randn(2,1); %starting point for fminsearch which is used for fitting

%extracting parameters for average expression
[b_average, i_average, d_average, tau_on_average, tau_off_average, maxl_average] = extract_cell(average, time, frames_base, frames_on, frames_off, frames_tail, x0, plotting_average, 1, 0);

%extracting paramters for single cells
[b, i, d, tau_on, tau_off, leakiness, maxl,  lmaxt] = extract(means, time, frames_base, frames_on, frames_off, frames_tail, x0);

%Note that the values used here are in arbitrary units obtained on the
%microscope 
%To convert them to promoter activity unit they should be divided by the maxGAL1 value
maxGAL1 = 5830.5; %maximal level of GAL1 expression

b_average_scaled = b_average/maxGAL1;
i_average_scaled = i_average/maxGAL1;
maxl_average = maxl_average/maxGAL1;

b_scaled = b/maxGAL1;
i_scaled = i/maxGAL1;
leakiness_scalled = leakiness/maxGAL1;
maxl_scalled = maxl/maxGAL1;
lmaxt_scalled = lmaxt/maxGAL1;
