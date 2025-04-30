Travelling wave project - to charactarize travelling waves in gamma (30-60Hz) from microelectrode recordings in NHP visual cortex. 
run 'runDisplayTWData.m'

The following programs should be downloaded and added to Matlab path 
1. Chronux: download from http://chronux.org/
2. CommonPrograms: https://github.com/supratimray/CommonPrograms
3. Circular statistics toolbox: https://github.com/circstat/circstat-matlab
4. Script for Watsons U2 test: https://github.com/pierremegevand/watsons_u2
** remove EEGLAB from path if added as it interferes with PCA in matlab (only for EEG data) **


Glossary of scripts
- load Data - load the unfiltered time series or reconstructed gamma data 
- getFilteredBurstsTW - filters the data and identifies gamma bursts based on hilbert method or sliding frequency as per Cohen MX. Fluctuations in oscillation frequency control spike timing and coordinate neural networks. Journal of Neuroscience. 2014 Jul 2;34(27):8988-98.
- getTWCircParams - run circular regression with burst detection for calculating TW parameters
- getHilbertBurst - modified script from gamma length project, to avoid filtering of data and feed reconstructed time series directly
- getWaveMetrics - a wrapper for circular linear regression (simple or cluster based) for grouping electrodes
- getWaveSegments - Segments/epochs the direction values to give a list of waves and their corrosponding time indices in the given trial
- circRegMod - calcualtes circular regression and give TW parameters
- getOverlappingWaves - identifies overlapping waves between slow and fast gamma in a given trial. 
 
Image Plots
- makePolarPLot - generate polar plots
- showRFPositionsSelectedElectrodes
- plotTFMT - plots the TF plots (delta or uncorrected) for a single channel, single/multiple trials given an orientation/SF combination
- displayTWData - GUI to view travelling wave parameters

Project Specific scripts (Dual Gamma Wave Project)
- generateFig1-5 - generate respective figures  
- runDisplayTWData - Gives an example of the entire pipeline

%%
Options in GUI
%% Data panel %%
Elecs - Specify 1-3 electrodes that are identified on the channel layout and RF plots. This determines the number
        of individual channel plots and time frequency plots. 
%% Axis Ranges 1 panel %%        
FreqLims - sets the y axis for the time freq plots
TimeLims - sets the x axis for all time Freq and line plots
cLims - sets colorLims for timefreq plots

%% Axis Ranges 2 panel %%
yRange (25-35Hz) - sets y axis ranges for all line plots corrosponding to slow gamma (on the left)
yRange (40-60Hz) - sets y axis ranges for all line plots corrosponding to fast gamma (on the right)

%% Travelling wave panel %% 
Electrode Fraction - '0.6' (default)
                   Once gamma bursts are detected, bursts across all electrodes are summed and only those time points
                   which have bursts across 60% (default) of electrodes are taken forward for wave detection. 
electrodeChoices - 'selected'/'all' 
                'selected' will use electrodes that show a gamma burst (detected by hilbert method or freequency sliding) and ignores all other electrodes.
                'all' uses all electrodes regardless of frequency content at the given time point.
 
wave Detection - 'Simple circLin regression'/'Cluster circLin regression'
               Simple circLin regression - Circular linear regression is done between the phase and location values at all selected electrodes
               Cluster circLin regression - Circular linear regression is done between the phase and location values only between clusters of electrodes. 
                                            Each cluster is considered to the nearest electrodes (within 2*interelectrode distance) of a given electrode.
                                            This method gives a unique direction for each electrode. Wave segmentation option is to be set to 'wave strength'. 

Wave wobble - '5' (default - in degree)
              This gives the limit of the variation in wave direction that is then used for wave segmentation   


waveSegOptions - 'Simple Segmentation'/'Point-point wobble'/'Full segment wobble'/'Wave strength'
                  Simple segmentation - waves are considered whererever significant PGD values (>0) are observed. Waves <25ms are ignored.
                  Point-Point Wobble - Based on the 'wave wobble', successive time points with significant PGD(>0) are said to belong to the 
                                       same wave only the wave direction varies by less than the 'wave wobble'.
                  Full segment wobble - A time segment is considered a wave if the variation between any two time points within
                                        the wave is less than the 'wave wobble'.
                  Wave strength - Determines the presence or absence of a wave at a given time point based on both the PGD and direction
                                  Details specified in Das, A., Zabeh, E., Jacobs, J. (2023). How to Detect and Analyze Traveling Waves in Human Intracranial EEG Oscillations?. In: Axmacher, N. (eds) Intracranial EEG. Studies in Neuroscience, Psychology and Behavioral Economics. Springer, Cham. https://doi.org/10.1007/978-3-031-20910-9_30
                                        
%% Plot Phase Propagation Panel %%
Time Range - set time range over which phase propagation is to be visualized
