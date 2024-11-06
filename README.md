Travelling wave project - to charactarize travelling waves in gamma (30-60Hz) from microelectrode recordings in NPH visual cortex. 
run 'runDisplayTWData.m'

The following programs should be downloaded and added to Matlab path 
1. Chronux: download from http://chronux.org/
2. CommonPrograms: https://github.com/supratimray/CommonPrograms
3. Circular statistics toolbox: https://github.com/circstat/circstat-matlab
4. Script for Watsons U2 test: https://github.com/pierremegevand/watsons_u2

Glossary of scripts
- runDisplayTWData - Gives an example of the entire pipeline
- displayTWData - GUI to view travelling wave parameters
- getTWCircParams - run circular regression with burst detection for calculating TW parameters
- getHilbertBurst - modified script from gamma length project, to avoid filtering of data and feed reconstructed time series directly
- load Data - load the unfiltered time series or reconstructed gamma data 
- circRegMod - calcualtes circular regression and give TW parameters
- getPolarPlor, getPolarPlotDos, getPolarPlorMean2 - generate polar plots
- findOverlappingWaves - gives waves in Slow and fast gamma that overlap by a specified threshold in all trials. Specify the type of segementation to be done. 
- getWaveSegments -  segments time series into waves. Does so by three methods, 1. Only consideres breaks in PGD to signify an individual wave 2. takes a relaxed 10degree deviation across each time point of a wave, 3. takes a stringent 5 degree deviation across the whole wave
- plotTFMT - plots the TF plots (delta or uncorrected) for a single channel, single/multiple trials given an orientation/SF combination
- generateFig1-5 - generate respective figures  

