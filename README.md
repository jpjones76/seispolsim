# SeisPolSim

## Introduction
This is the online repository for programs based on `[1]`. I will maintain and improve these programs as time allows. At this point, code is stable but rigid. Programs based on later iterations of `[2]` and subsequent research will be added as the project develops.

# Polarization Similarity
The following sections describe how to use the polarization similarity routines. Additional information can be found in the help text of various files.

## seispol
The most basic program for computing polarization similarity. Workflow is described in `[1]`. An example run might look like this:

1. Preprocess first. Rudimentary preprocessing (e.g. de-mean, instrument response) can be done with scripts in Misc, if needed.
2. Load data into matrix **X** in column vectors.
  * One column per channel.
  * Arrange data in order Z,N,E.
  * Treat +Z as down.
  * Expected column order is e.g. `[z1 n1 e1 ... zK nK eK]`. 
3. `L = 100; La = 50;` 
  * Initialize some useful variables. 
  * We'll use L for the length of each histogram in samples.
  * We'll use La for the number of samples to advance between histogram calculations.
4. `[P, W, H1] = seispol(X, 'dp', 0.01, 'h1', 1);`
  * `'h1', 1`: Compute a single histograms for each channel. Useful if input data contain a single phase or sensor alignment is a concern.
  * `'dp', 0.01`: Normalized bin spacing; each histogram has 1/dp bins.
5. `[~, ~, H2, D, T] = seispol(X, 'L', L, 'La', La, 'dp', 0.01, 'ht', 1);`
  * `'ht', 1`: Compute time-slice histograms. 
  * `'L', L`: Histogram windows will each be L samples long. 
  * `'La', La`: Advance by La samples between histograms.
  * `'dp', 0.01`: Normalized bin spacing; each histogram has 1/dp bins.
6. `S = adaptivesim(X, D, L, La);`
  * Convert distances to adaptive similarities. 

Example Script: Demos/eq_banff_polsim.m

## polsim_modwt
This is the most basic time-frequency analysis routine: Transform to the discrete wavelet domain via undecimated DWT (aka "MODWT"), create detail coefficients, and compute polarization of those coefficients. Largely for illustrative purposes. Future work will use either the undecimated DWPT, as in [2], or the STFT, as appropriate.

Example Script: Demos/eq_banff_modwt.m

## Visualization
A few scripts have been included in Graphics. These are meant to be better-organized versions of material used in [1-2].

### psiplot
Create intensity plots of polarization similarity. Two plots are generated.

1. Adjacent-time simlarity: For attribute p, average similarity in each time slice t between `H{p}(t,k) H{p}(t,k1~=k)`.
2. Cross-station similarity: For attribute p, similarity at each station k between `H{p}(t,k), H{p}(t-1,k)`.

Example use: `h = psiplot(X,S,T,fs,L);`

###pstfip
As psiplot, for the output of polsim_modwt. Four sets of intensity plots are created, each with octave scaling on the frequency axis.

1. Cross-station similarity: For attribute p, similarity at each station k between `H{p}(t,k), H{p}(t-1,k)`, in each set of MODWT coefficients Wj. Each band has a ragged, overloaded look because they contain K total sets of intensities, stacked vertically. These plots are meant for diagnostic purposes.
2. Adjacent-time simlarity: For attribute p, max. MODWT level J, average similarity in each time slice t between `H{p}(t,k) H{p}(t,k1~=k)` at each level j of the wavelet transform. As above, the look is ragged because each band has K sets of intensities.
3. Averaged cross-station similarity.
4. Averaged adjacent-time similarity for all stations. As (1), but averaged across all stations considered at each scale j.

###Tips
* Use common sense when determining window lengths. For example, 2s histograms are too short for teleseisms, too long for very small transients, but quite good for regional earthquakes. 

## Data Availability
Sample data used in our publications may be made available upon request. We do not have permission to distribute it freely. Please follow the guidelines below.

### Hoadley Borehole Data
All requests must go through [Prof. David W. Eaton](http://www.ucalgary.ca/eatond/). Please contact J. Jones if you have trouble reaching Prof. Eaton.

### Yukon/BC Data
Data from all permanent CNSN stations are publicly available from [Earthquakes Canada AutoDRM](http://www.earthquakescanada.nrcan.gc.ca/stndon/AutoDRM/index-eng.php). 

### Paper Scripts and Plotting Utils
Scripts to (re)generate the figures in `[1-2]` are available from J. Jones by request.

## Functionality
The initial set of programs should be usable (without modification) in Matlab/Octave; download all and run install.m to compile the .mex files. 

### Note to Apple Users
This project currently has no Apple support for time-frequency analysis. 

### Issues with Plotting
Graphics-related routines may have backwards compatibility issues with Matlab < r2014b. 

## References

1. Jones, J.P., Eaton, D., and Caffagni, E., 2015, "Quantifying the similarity of seismic polarizations", *GJI, in press*.
2. Jones, J.P., Eaton, D., and Caffagni, E., 2015, "Quantifying similarity in seismic polarizations", *AGU Fall Meeting*, Paper ID S21B-2694.
