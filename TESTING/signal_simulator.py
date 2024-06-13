import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from scipy.stats import linregress
from fractions import Fraction
import _LCurveSimulator as cppLCurveSim # import C++ script to simulate signals
PATH0 = os.getcwd()
sampling_files = [file for file in os.listdir() if file.endswith(".dat")]

## LATEX FONTS:
if True:
  plt.rcParams['text.latex.preamble']= r"\usepackage{lmodern}"
  params = {'text.usetex' : True,
				'font.size' : 21,
				'font.family' : 'lmodern'}
  plt.rcParams.update(params)



SAMPLINGS = ['1_3C454_Medicina_X', '2_0716_fgamma_6cm', '3_0716_urumqi_200607', '4_0212_gbi', '5_0133_gbi', '6_0235_gbi', '7_even']
DIVIDE_BETA = False # True to generate signals with fixed beta values for the Power Spectrum and divide the signals consequently; if False, will generate signals with an expanded sample of values of beta
if DIVIDE_BETA:
	NSIMULATIONS = 1111 ### number of simulated signals
	BETAS = [-1.5,-2.,-2.5] # fixed amount of betas
else:
	NSIMULATIONS = 333 ### number of simulated signals
	Nbetas = 12 ### will get simulations for Nbetas-1 slopes of the PSD
	BETAS = np.unique(np.round(np.concatenate((np.linspace(-2.5,-2.,int(Nbetas/2)), np.linspace(-2.,-1.5,int(Nbetas/2)))), 2))
	### If we need to randomly generate beta values:
	#random_betas = np.round(np.random.uniform(-2.5, -1.5, Nbetas-3),2)
	#BETAS = np.concatenate([random_betas, [-1.5, -2., -2.5]]) # ensure a sample of interesting beta values, p.e. -2: typical red noise signal
	### Write in a .txt file for later uses
	#with open("beta_values.txt", "w") as file:
	#	file.write(", ".join(map(str, BETAS)))
	### To read the file:
	### with open("beta_values.txt", "r") as file:
	###	BETAS = list(map(float, file.read().split(", ")))

NOISE_LEVEL = [3.]

def LCurveSimulator(times,beta,NOISElvl=[3.],Nwave=50,fname='simulated.dat'):
	'''
		Function to simulate a Light Curve
		Parameters:
		- times: array with the Light Curve times
		- beta: power law spectrum factor
		- noiselvl: percentage of the simulated signal to add as noise
		- Nwave: increase "Nwave/2" times the "Ntimes" of the original Signal
		- fname: name of the output ascii file with the simulated LCurve
	'''
	Ntimes = len(times)
	### Step 0: get array with frequency values, given by signal in Fourier domain
	TOT_TIME = np.max(times) - np.min(times)
	freqs = np.array([n/(Nwave*TOT_TIME) for n in np.arange(1,int(Nwave*len(times)/2.))])
	Nfreqs = len(freqs)
	### Prepare times and freqs for the C++ function
	times_cpp = np.require(np.copy(times),requirements=['C','A'])
	freqs_cpp = np.require(np.copy(freqs),requirements=['C','A'])
	### Initialize output variables for the C++ function
	FTreal = np.require(np.zeros(Nfreqs,dtype=np.float64),requirements=['C','A'])
	FTimag = np.require(np.zeros(Nfreqs,dtype=np.float64),requirements=['C','A'])
	signal = np.require(np.zeros(Ntimes,dtype=np.float64),requirements=['C','A'])
	
	cppLCurveSim.SimSignalcpp(times_cpp,freqs_cpp,Ntimes,Nfreqs,beta,FTreal,FTimag,signal)
	FT = FTreal * FTreal + FTimag * FTimag ### it is like this because FTreal and FTimag are random gaussian numbers of amplitude sqrt(0.5*S_freq)
	
	fig = pl.figure(figsize=(10,7))
	sub1 = fig.add_subplot(111)
	fig.subplots_adjust(wspace=0.01,hspace=0.01,right=0.98,left=0.125)
	fig.suptitle(fname[:-4] + ' beta=%.2f'%beta,fontsize=21)
	sub1.loglog(freqs,FT,'.k', label='Spectrum')
	slope, intercept, _, _, _ = linregress(np.log(freqs), np.log(np.abs(FT)))
	sub1.loglog(freqs, np.exp(intercept) * freqs**slope, label=f'Slope: {slope:.3f}, logT0: {intercept:.3f}') ### linear fitting to log(FT) vs log(freqs)
	sub1.set_xlabel(r'Frequency (h$^{-1}$)')
	sub1.set_ylabel('Amplitude')
	pl.legend()
	pl.savefig(fname[:-4]+'Fourier.png')
	pl.close()
	
	### Save and plot the signal in time domain, i.e. after FFT (Fast Fourier Transform)
	for noiselvl in NOISElvl:
		noise = np.random.normal(0, noiselvl/100.*np.std(signal), len(signal))
		signal_contaminated = signal+noise
		with open(fname[:-4]+"_noise%i.dat"%noiselvl, "w+") as fpOut:
			for ti, amp, err in zip(times,signal_contaminated,noise):
				fpOut.write("{:.7e}   {:.7e}   {:.7e}\n".format(ti,amp,err))
		fig = pl.figure(figsize=(10,7))
		sub1 = fig.add_subplot(111)
		fig.subplots_adjust(wspace=0.01,hspace=0.01,right=0.98,left=0.125)
		fig.suptitle(fname[:-4] + ' noise=%i beta=%.2f'%(noiselvl,beta),fontsize=21)
		sub1.plot(times,signal_contaminated,'.k')
		sub1.set_xlabel('Times (h)')
		sub1.set_ylabel('Amplitude')
		pl.savefig(fname[:-4]+'_noise%i.png'%noiselvl)
		pl.close()
	del Ntimes,TOT_TIME,freqs,Nfreqs,times_cpp,freqs_cpp,FTreal,FTimag,signal,FT,noise,signal_contaminated,slope,intercept

import time
init_time = time.time()

for sampling in SAMPLINGS:
	dataname = 'nsampling%s.dat'%sampling
	### create directory to store simulated signals for a given sampling
	PATH_sampling = os.path.join(PATH0,sampling)
	if not os.path.exists(PATH_sampling):
		os.makedirs(PATH_sampling)
	### load times or skip is missing
	if dataname in sampling_files:
		data = np.loadtxt(dataname)
		Time = data[:,0]
		del data
	else:
		print("Sampling not found in the current directory\n")
		continue
	os.chdir(PATH_sampling)
	
	for beta in BETAS:
		print('\nSimulating Signal for sampling%s with beta=%.2f\n'%(sampling,beta))
		### run simulations
		for nsim in range(NSIMULATIONS):
			if nsim%10==0:
				sys.stdout.write('\rDoing simulation %i of %i\n'%(nsim,NSIMULATIONS))
			fOutname = 'nsampling%s_beta%.2f_sim%i.dat'%(sampling,beta,nsim)
			LCurveSimulator(Time,beta,NOISE_LEVEL,50,fOutname)
		os.chdir(PATH_sampling)
	os.chdir(PATH0)

end_time = time.time()
totat_time_duration = end_time - init_time
print("Total exec. time: {:.4f} seconds.".format(totat_time_duration))

