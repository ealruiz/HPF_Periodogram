import glob, os, gc, re
import numpy as np
import matplotlib.pyplot as plt
PATH0 = os.getcwd()

## LATEX FONTS:
if True:
  plt.rcParams['text.latex.preamble']= r"\usepackage{lmodern}"
  params = {'text.usetex' : True,
				'font.size' : 21,
				'font.family' : 'lmodern'
				}
  plt.rcParams.update(params)

SAMPLINGS = ['1_3C454_Medicina_X', '2_0716_fgamma_6cm', '3_0716_urumqi_200607', '4_0212_gbi', '5_0133_gbi', '6_0235_gbi', '7_even']
Nsamplings = len(SAMPLINGS)
DIVIDE_BETA = True # True to generate signals with fixed beta values for the Power Spectrum and divide the signals consequently; if False, will generate signals with an expanded sample of values of beta
if DIVIDE_BETA:
	NSIMULATIONS = 1111 ### number of simulated signals
	BETAS = [-1.5,-2.,-2.5] # fixed amount of betas
	Nbetas = len(BETAS)
else:
	NSIMULATIONS = 33 ### number of simulated signals
	Nbetas = 7
	BETAS = np.round(np.concatenate((np.linspace(-2.5,-2.,int(Nbetas/2)), np.linspace(-2.,-1.5,int(Nbetas/2)))), 2)

NOISE_LEVEL = [3.]
Nnoise = len(NOISE_LEVEL)

METHODS = ['Dynamic','BottomUp', 'Kernlin'] # 'Binary', 'Kernrbf',

### Parameters of the simulated periodic signal
X0_FACT = [0.3,1.,3.] # the amplitude of the periodic signal will be X0_FACT times the standard deviation of the random simulated signal
PERIOD_LEN = [0.01, 0.05, 0.1, 0.25, 0.4] # the period of the periodic signal (freq0=1/period) will be PERIOD_LEN times the total lenght (i.e. the duration) of the simulated signal
### Parameters of the simulated periodic signal, inlcuding 0 (original signal without periodicity)
ALL_X0_FACT = [0.] + X0_FACT
NX0 = len(ALL_X0_FACT)
ALL_PERIOD_LEN = [0.] + PERIOD_LEN
NPERIOD = len(ALL_PERIOD_LEN)
VERSIONS = ['original','detrended']

P_false_alarm = 0.01 # False alarm probability for periodicity detections

def GetPeriodicityDetections(fname='normalised_periodogram.dat',FreqPeriod=0.):
	data = np.loadtxt(fname,skiprows=1)
	freqs = np.array([ti for ti in data[:,0]])
	Nfreqs = len(freqs)
	norm_pgram = np.array([si for si in data[:,1]])
	z0 = np.log(Nfreqs/P_false_alarm)
	Ndetections, PeriodicityFound, Periodicity_MaxPeriodogram = 0, 0, 0
	for indx,npgram in enumerate(norm_pgram):
		if z0 < npgram:
			Ndetections += 1
			if indx>0 and indx<(Nfreqs-1) and (freqs[indx-1]<=FreqPeriod<=freqs[indx+1]):
				PeriodicityFound = 1
				if PeriodicityFound==1 and npgram==np.max(norm_pgram):
					Periodicity_MaxPeriodogram = 1
	return Ndetections, PeriodicityFound, Periodicity_MaxPeriodogram

######## STEP 1: locate possible periodicities

with open("Periodicity_Detections.txt", "w") as file:
		file.write("Sampling\tVersion\tbeta\tNoise\tX0\tPeriod\tMethod\tDetections\tPeriodicity Found?\tPeriodicity at Max. NormPeriodogram?\n")

for sampling in SAMPLINGS:
	datname = 'nsampling%s.dat'%sampling
	simsignal = np.loadtxt(datname)
	TOT_TIME = np.max(simsignal[:,0]) - np.min(simsignal[:,0])
	PATH_sampling = os.path.join(PATH0,sampling)
	for version in VERSIONS:
		if version == 'original':
			PATH_PERIOD = os.path.join(PATH_sampling,'PERIODOGRAMS')
		else:
			PATH_PERIOD = os.path.join(PATH_sampling,'DETRENDED_PERIODOGRAMS')
		for noiselvl in NOISE_LEVEL:
			for beta in BETAS:
				for method in METHODS:
					print('\nSearch periodicity in %s %s noise%s beta%s %s\n'%(sampling,version,noiselvl,beta,method))
					os.chdir(PATH_PERIOD)
					FreqPeriod = 0.
					Ndetections, PeriodicityFound, Periodicity_MaxPeriodogram = 0, 0, 0
					for nsim in range(NSIMULATIONS):
						if version == 'original':
							fname = 'nsampling%s_beta%.2f_sim%i_noise%i_periodogram_normalised_%s.dat'%(sampling,beta,nsim,noiselvl,method)
						else:
							fname = 'nsampling%s_beta%.2f_sim%i_noise%i_periodogram_detrended_normalised_%s.dat'%(sampling,beta,nsim,noiselvl,method)
						n_detections, periodicity_found, max_pgram_period = GetPeriodicityDetections(fname,FreqPeriod)
						Ndetections += n_detections
						PeriodicityFound += periodicity_found
						Periodicity_MaxPeriodogram += max_pgram_period
					os.chdir(PATH0)
					with open("Periodicity_Detections.txt", "a") as file:
						file.write("%s\t%s\t%s\t%s\t%.1f\t%.2f\t%s\t%i\t%i\t%i\n"%(sampling,version,beta,noiselvl,0.,0.,method,Ndetections,PeriodicityFound, Periodicity_MaxPeriodogram))
					for amp_fact in X0_FACT:
						PATH_X0 = os.path.join(PATH_sampling,'PSignal_%.1f'%amp_fact)
						if version == 'detrended':
							PATH_X0_PERIODOGRAM = os.path.join(PATH_X0,'Detrended_NormPeriodogram')
						else:
							PATH_X0_PERIODOGRAM = os.path.join(PATH_X0,'NormPeriodogram')
						for period in PERIOD_LEN:
							os.chdir(PATH_X0_PERIODOGRAM)
							FreqPeriod = 1/(TOT_TIME*period)
							Ndetections, PeriodicityFound, Periodicity_MaxPeriodogram = 0, 0, 0
							for nsim in range(NSIMULATIONS):
								if version == 'detrended':
									fname = 'nsampling%s_beta%.2f_sim%i_noise%i_Xo%.1f_T%.2f_periodogram_detrended_%s_normalised_%s.dat'%(sampling,beta,nsim,noiselvl,amp_fact,period,method,method)
								else:
									fname = 'nsampling%s_beta%.2f_sim%i_noise%i_Xo%.1f_T%.2f_periodogram_%s_normalised_%s.dat'%(sampling,beta,nsim,noiselvl,amp_fact,period,method,method)
								n_detections, periodicity_found, max_pgram_period = GetPeriodicityDetections(fname,FreqPeriod)
								Ndetections += n_detections
								PeriodicityFound += periodicity_found
								Periodicity_MaxPeriodogram += max_pgram_period
							os.chdir(PATH0)
							with open("Periodicity_Detections.txt", "a") as file:
								file.write("%s\t%s\t%s\t%s\t%.1f\t%.2f\t%s\t%i\t%i\t%i\n"%(sampling,version,beta,noiselvl,amp_fact,period,method,Ndetections,PeriodicityFound, Periodicity_MaxPeriodogram))
