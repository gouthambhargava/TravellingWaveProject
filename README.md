Travelling wave project - to charactarize travelling waves in gamma (30-60Hz) from microelectrode recordings in NPH visual cortex. 
run 'runDisplayTWData.m'

The following programs should be downloaded and added to Matlab path 
1. Chronux: download from http://chronux.org/
2. ComonPrograms: https://github.com/supratimray/CommonPrograms
3. Gamma Length Project repository: https://github.com/supratimray/GammaLengthProjectCodes)
4. Supporting matlab files for MP: http://www.fuw.edu.pl/~durka/software/mp/.
5. Cricular statistics toolbox: https://github.com/circstat/circstat-matlab

Glossary of scripts
- runDisplayData - Gives an example of the entire pipeline
- displayTWData - GUI to view travelling wave parameters
- getTWCircParams - run circular regression with burst detection for calculating TW parameters
- RunMP - run matching pursuit (stochastic algo) to extract gabor data and to save reconstructed time series
- getHilbertBurst - modified script from gamma length project, to avoid filtering of data and feed reconstructed time series directly
- load Data - load the unfiltered time series or reconstructed gamma data 
- circRegMod - calcualtes circular regression and give TW parameters
- gammaBurstMPParams - extracts data following MP burst detection for plotting

