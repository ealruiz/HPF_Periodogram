# New periodogram implementation

This project includes the scripts needed to compute the periodogram of a signal using a new implementation proposed by Ezequiel Albentosa Ruiz and Nicola Marchili. Together with the main scripts, a pipeline is given for testing the performance of this new implementation compared with the classic periodogram proposed by Jeffrey D. Scargle (1982).

The contents of the repository are explained in the following sections:

## Main scripts:
- detrendedPeriodogram.cpp: C++ module to compute the new implementation of the periodogram. This script has to be imported in python. This script requires the following arguments from python:
  - Time, LCurve, LCurve_err: basic data of the signal: time series, and signal amplitudes and measurement errors.
  - PFreq, Periodogram: frequencies of the signal in Fourier domain, at which the periodogram will be evaluated, and array, with the same lenght as PFreq, to be filled with the periodogram values.
  - Ndata, Nfreqs: number of data in the time series and number of frequencies.
  - LogScale: default: 0; set to 1 to compute the log. scale of the signal; useful for signals with large spikes (p.e. optical lcurves or pulsar detections).
  - save_detrend_signal: default: 0; set to 1 to save the detrended signals for each tbin (maybe useful for testing the working of the periodogram).
- detrendedPeriodogram_setup.py: execute this python script to setup C++ module "detrendedPeriodogram.cpp".

## Testing performance
Grouped in the "TESTING" folder, you will find the following scripts needed for testing the performance of the new periodogram implementation and compare with the classic periodogram:
- "*.dat" files with the time series for different samplings. The amplitudes correspond to a sinusoidal signal, but are not used.
- LCurveSimulator.cpp: C++ module to simulate signals. This script has to be imported in python. This script requires the following arguments from python:
  - times, freqs: time series of the signals, and frequencies in Fourier domain.
  - beta: Power Spectrum Distribution (PSD) of the signal in Fourier domain.
  - Noise_Level: percentage of the simulated signal std to add as noise
  - Ndata, Nfreqs: number of times in the time series and number of frequencies.
  - FTreal,FTimag: FT (real and imaginary part) of the simulated signal in Fourier domain.
  - signal: simulated signal amplitudes.
- LCurveSimulator_setup.py: execute this python script to setup C++ module "LCurveSimulator.cpp".
- ComputePeriodogram.cpp: C++ module to compute the classic periodogram. This script has to be imported in python. This script requires the following arguments from python:
  - Time, LCurve, LCurve_err: basic data of the signal: time series, and signal amplitudes and measurement errors.
  - PFreq, Periodogram: frequencies of the signal in Fourier domain, at which the periodogram will be evaluated, and array, with the same lenght as PFreq, to be filled with the periodogram values.
  - Ndata, Nfreqs: number of data in the time series and number of frequencies.
  - LogScale: default: 0; set to 1 to compute the log. scale of the signal; useful for signals with large spikes (p.e. optical lcurves or pulsar detections).
- ComputePeriodogram_setup.py: execute this python script to setup C++ module "ComputePeriodogram.cpp".

The following scripts will use the C++ modules and some parameters to simulate signals, compute the periodogram, and proceed with the testing: (they have to be executed sequently):
1. signal_simulator.py: simulates the signals using the "LCurveSimulator.cpp" module. Will save the simulated signals in a separate folder for each sampling.
2. periodogram_detrended.py: python script to compute the new implementation of the periodogram using the "detrendedPeriodogram.cpp" module. Will save the periodograms in a separate folder for each sampling.
3. periodogram.py: python script to compute the classic periodogram using the "ComputePeriodogram.cpp" module. Will save the periodograms in a separate folder for each sampling.
4. Periodogram_Fit.py: script to get PSD from the Classic and New Periodogram of the simulated signals, and compare with the expected (simulated) PSD value. This version should be used to focus on testing PSD estimates. Will save all results in a "periodogram_slopes.txt" file, with the following columns: Sampling, beta (simulated PSD), Nsim (number to identify each simulation), Noise (noise level with respect the simulated signal std), Method (method used to get noise level), Slope (PSD estimate from the periodogram), Err. (discrepancy of the PSD estimate with respect the simulated PSD), and Version ("original" for the classic periodogram, "detrended" for the new implementation).
5. Periodogram_Periodicity_Fit_and_Normalise.py: script to add periodic signal to the simulated signals, fit (and normalise) the periodograms (both classic and new implementation) of the original signal first, and then of the signals with periodicity. This script should be used to focus on testing periodicity detections. Will save the normalised periodograms of each simulation (needed to search for periodicity detections), all PSD results in a "periodogram_slopes.txt" file (with same columns as the previous step), and an additional "periodic_signal_periodogram_slopes.txt" file with all the PSD results , with the following columns (including two additons to identify the periodic signal): Sampling, beta, Nsim, Noise, X0 (to identiry the amplitude of the periodic signal), Period (to identify the period of the periodic signal), Method, Slope, Err., Version.
6. Periodicity_statistics_detection.py: script to locate possible periodicities from the normalised periodograms obtained with the "Periodogram_Periodicity_Fit_and_Normalise.py" script. Will save the results in a "Periodicity_Detections.txt"" file, with the following columns (for each group of Nsim): Sampling, Version, beta, Noise, X0, Period, Method, Detections (number of frequencies at which we can claim detection using the periodogram), Periodicity Found? (number of simulations at which we find periodicity at the expected frequency, i.e. at freq0 = 1/period), Periodicity at Max. NormPeriodogram? (number of simulations at which the periodicity found at the expected frequency corresponds to the maximum value of the periodogram).
- These scripts require the following parameters:
  - SAMPLINGS: array with the names of the samplings used to simulate the signals
  - NSIMULATIONS: number of simulations to generate.
  - BETAS: array with the PSD of the signals in Fourier domain. Will simulate NSIMULATIONS per beta value.
  - NOISE_LEVEL: array with the percentage of the simulated signal std to add as noise. Only 3% level by default (i.e. NOISE_LEVEL=[3.]), but more values can be added for testing. Will simulate NSIMULATIONS per beta value.
- Other parameters, only needed for some of the scripts, are:
  - PLOT_PERIODOGRAM: default: False (only needed for the periodogram scripts); set True to plot the periodogram (will increase execution time)
  - LogScale: default: 0; set to 1 to compute the log. scale of the signal (only needed for the periodogram scripts)
  - save_detrend_signal: default: 0 (only needed for the "periodogram_detrended.py" script); set to 1 to save the detrended signals for each tbin (will increase execution time)
  - METHODS: array with methods to get noise level automatically, by identifying the break in the data trend (from PSD to flat). Only needed for the "Periodogram_Fit.py" script. Recommended: ['BottomUp', 'Binary', 'Dynamic', 'Kernrbf', 'Kernlin']. See documentation of the python "ruptures" library for more info.
  - X0_FACT: array with the amplitude factor of the periodic signal; the amplitude of the periodic signal will be X0_FACT times the standard deviation of the random simulated signal
  - PERIOD_LEN: array with the period length factor; the period of the periodic signal (freq0=1/period) will be PERIOD_LEN times the total lenght (i.e. the duration) of the simulated signal. Only needed for the "Periodogram_Periodicity_Fit_and_Normalise.py" script.

