import glob, os, gc, sys
import numpy as np
import math
import matplotlib.pyplot as plt
import pylab as pl
from fractions import Fraction
import _ComputePeriodogram as Periodogramcpp # import C++ script to compute Classic Periodogram
PATH0 = os.getcwd()

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

NOISE_LEVEL = [3.]

PLOT_PERIODOGRAM = False ### True to plot the periodogram
LogScale = 0 ### set to 1 to compute the log. scale of the signal; useful for signals with large spikes (p.e. optical lcurves)

def GetPeriodogram(fname='periodogram.dat',plotfig=False):
	data = np.loadtxt(fname)
	times = data[:,0]
	Ndata = len(times)
	TOT_TIME = np.max(times) - np.min(times)
	freqs = np.array([2.*np.pi*n/TOT_TIME for n in np.arange(1,int((len(times)+1)/2))])
	Nfreqs = len(freqs)
	### require copy needed of selected data column for C++ processing
	Time = np.require(np.copy(times),requirements=['C','A'])
	PFreq = np.require(np.copy(freqs),requirements=['C','A'])
	if LogScale == 1:
		LCurve = np.require(np.copy(np.log10(data[:,1])),requirements=['C','A'])
	else:
		LCurve = np.require(np.copy(data[:,1]),requirements=['C','A'])
	os.chdir(PATH_PERIOD)
	Periodogram = np.require(np.zeros(Nfreqs,dtype=np.float64),requirements=['C','A'])
	### Get periodogram
	Periodogramcpp.ComputePeriodogram(Time,LCurve, PFreq,Periodogram, Ndata,Nfreqs)
	### Save and  Plot the periodogram
	Outdata = np.column_stack((freqs/(2.*np.pi), Periodogram))
	np.savetxt(fname[:-4]+'_periodogram.dat', Outdata, delimiter='\t', header='#\t Freq (h^-1) \t Periodogram', comments='')
	if plotfig:
		fig = pl.figure(figsize=(10,7))
		sub1 = fig.add_subplot(111)
		fig.subplots_adjust(wspace=0.01,hspace=0.01,right=0.98,left=0.125)
		fig.suptitle(fname[:-4],fontsize=21)
		sub1.loglog(np.array(PFreq)/(2.*np.pi),Periodogram,'.k')
		sub1.set_xlabel(r'Frequency (h$^{-1}$)')
		sub1.set_ylabel('Periodogram')
		pl.savefig(fname[:-4]+'_periodogram.png')
		pl.close()
	del data,times,Ndata,TOT_TIME,freqs,Nfreqs,Time,PFreq,LCurve,Periodogram,Outdata



import time
init_time = time.time()


for sampling in SAMPLINGS:
	PATH_sampling = os.path.join(PATH0,sampling)
	os.chdir(PATH_sampling)
	for beta in BETAS:
		PATH_PERIOD = os.path.join(PATH_sampling,'PERIODOGRAMS')
		if not os.path.exists(PATH_PERIOD):
			os.makedirs(PATH_PERIOD)	
		for nsim in range(NSIMULATIONS):
			for noiselvl in NOISE_LEVEL:
				print('\nComputing periodogram for sampling%s beta=%s SIM%i noiselvl=%s\n'%(sampling,beta,nsim,noiselvl))
				fname = 'nsampling%s_beta%.2f_sim%i_noise%i.dat'%(sampling,beta,nsim,noiselvl)
				GetPeriodogram(fname,PLOT_PERIODOGRAM)
				os.chdir(PATH_sampling)


end_time = time.time()
totat_time_duration = end_time - init_time
print("Total exec. time: {:.4f} seconds.".format(totat_time_duration))

