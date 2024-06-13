import glob, os, gc, sys
import numpy as np
import math
import matplotlib.pyplot as plt
import pylab as pl
from fractions import Fraction
import _detrendedPeriodogram as detrendPeriodogramcpp # import C++ script to compute Detrended Periodogram
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
save_detrend_signal = 0 ### set to 1 to save the detrended signals for each tbin


def GetDetrendPeriodogram(fname='periodogram.dat',plotfig=False):
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
	#LCurve_err = np.require(np.copy(data[:,2]),requirements=['C','A'])
	LCurve_err = np.zeros(len(data[:,1])) ### careful: no error used!
	os.chdir(PATH_PERIOD)
	Periodogram = np.require(np.zeros(Nfreqs,dtype=np.float64),requirements=['C','A'])
	### Get periodogram
	detrendPeriodogramcpp.detrendedPeriodogram(Time,LCurve,LCurve_err, PFreq,Periodogram, Ndata,Nfreqs, LogScale, save_detrend_signal)
	### Sabe and Plot the periodogram
	Outdata = np.column_stack((freqs/(2.*np.pi), Periodogram))
	np.savetxt(fname[:-4]+'_periodogram_detrended.dat', Outdata, delimiter='\t', header='#\t Freq (h^-1) \t Periodogram', comments='')
	if plotfig:
		fig = pl.figure(figsize=(10,7))
		sub1 = fig.add_subplot(111)
		fig.subplots_adjust(wspace=0.01,hspace=0.01,right=0.98,left=0.125)
		fig.suptitle(fname[:-4],fontsize=21)
		sub1.loglog(np.array(PFreq)/(2.*np.pi),Periodogram,'.k')
		sub1.set_xlabel(r'Frequency (h$^{-1}$)')
		sub1.set_ylabel('Periodogram')
		pl.savefig(fname[:-4]+'_periodogram_detrended.png')
		pl.close()
	
	if save_detrend_signal==1: ### plot all intermediate interp. used for periodogram
		for indx in range(len(freqs)):
			fdatname = "detrended_signal_freq%i.dat"%indx
			if indx%111==0:
				sys.stdout.write("\rReading and ploting: %s\n"%fdatname)
			detrend_signal_data = np.loadtxt(fdatname)
			spline_signal = detrend_signal_data[:,1]
			detrended_signal = detrend_signal_data[:,2]
			plt.figure(figsize=(10,7))
			plt.subplots_adjust(wspace=0.01,hspace=0.01, right=0.99,left=0.125,top=0.975)
			plt.scatter(times, LCurve, label='Signal')
			plt.scatter(times, spline_signal, label='Spline signal')
			plt.scatter(times, detrended_signal, label='Detrend. signal')
			plt.legend()
			plt.savefig(fdatname[:-4]+".png")
			plt.close()
	del data,times,Ndata,TOT_TIME,freqs,Nfreqs,Time,PFreq,LCurve,LCurve_err,Periodogram,Outdata




import time
init_time = time.time()

for sampling in SAMPLINGS:
	PATH_sampling = os.path.join(PATH0,sampling)
	os.chdir(PATH_sampling)
	for beta in BETAS:
		PATH_PERIOD = os.path.join(PATH_sampling,'DETRENDED_PERIODOGRAMS')
		if not os.path.exists(PATH_PERIOD):
			os.makedirs(PATH_PERIOD)	
		for nsim in range(NSIMULATIONS):
			for noiselvl in NOISE_LEVEL:
				print('\nComputing detrended periodogram for sampling%s beta=%s SIM%i noiselvl=%s\n'%(sampling,beta,nsim,noiselvl))
				fname = 'nsampling%s_beta%.2f_sim%i_noise%i.dat'%(sampling,beta,nsim,noiselvl)
				GetDetrendPeriodogram(fname,PLOT_PERIODOGRAM)
				os.chdir(PATH_sampling)


end_time = time.time()
totat_time_duration = end_time - init_time
print("Total exec. time: {:.4f} seconds.".format(totat_time_duration))

