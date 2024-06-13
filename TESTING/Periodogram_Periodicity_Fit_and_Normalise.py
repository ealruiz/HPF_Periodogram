import glob, os, gc, re
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from scipy.optimize import curve_fit
from fractions import Fraction
import ruptures as rpt
import _detrendedPeriodogram as detrendPeriodogramcpp # import C++ script to compute Detrended Periodogram
import _ComputePeriodogram as Periodogramcpp # import C++ script to compute Classic Periodogram
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
DIVIDE_BETA = True # True to generate signals with fixed beta values for the Power Spectrum and divide the signals consequently; if False, will generate signals with an expanded sample of values of beta
if DIVIDE_BETA:
	NSIMULATIONS = 1111 ### number of simulated signals
	BETAS = [-1.5,-2.,-2.5] # fixed amount of betas
else:
	NSIMULATIONS = 33 ### number of simulated signals
	Nbetas = 7
	BETAS = np.round(np.concatenate((np.linspace(-2.5,-2.,int(Nbetas/2)), np.linspace(-2.,-1.5,int(Nbetas/2)))), 2)

NOISE_LEVEL = [3.]

METHODS = ['Dynamic','BottomUp', 'Kernrbf'] # 'Binary', 'Kernlin',
PLOT_FIT = False # True to plot "linear" fit of the Periodogram vs freq. in log-log scale

### Parameters of the simulated periodic signal to be injected
X0_FACT = [0.3,1.,3.] # the amplitude of the periodic signal will be X0_FACT times the standard deviation of the random simulated signal
PERIOD_LEN = [0.01, 0.05, 0.1, 0.25, 0.4] # the period of the periodic signal (freq0=1/period) will be PERIOD_LEN times the total lenght (i.e. the duration) of the simulated signal


# Define the functions for the fits
def linear_fit(x, a, b):
	return a * x + b


def getBreakIndx(data,method,MinData=11,NBreaks=1,PlotBreaks=False):
	"""
		Function with methods to search for the change points in a series of data points
		Parameters:
			- data: data points
			- method: method to get the breaks in the data (description at the end of the script)
			- MinData: minimum number of data points in the segment
			- NBreaks: number of change points in the data
			- PlotBreaks: set to True to plot the breaks in the data
		Return: 
			- result: array of indices where there is a change point in the data series
	"""
	method_algorithms = {
		 "Dynamic": rpt.Dynp(model="l2", min_size=MinData),
		 "Binary": rpt.Binseg(model="l2", min_size=MinData),
		 "BottomUp": rpt.BottomUp(model="l2", min_size=MinData),
		 "Kernlin": rpt.KernelCPD(kernel="linear", min_size=MinData),
		 "Kernrbf": rpt.KernelCPD(kernel="rbf", min_size=MinData),
		 "Kerncos": rpt.KernelCPD(kernel="cosine", min_size=MinData)
	}

	if method in method_algorithms:
		algo = method_algorithms[method]
		algo.fit(data)
		result = algo.predict(n_bkps=NBreaks) # number of change points in the data
	else:
		print('\nWARNING: method %s not in list of methods\n'%method)
		result = [1,len(data)-1]
	if PlotBreaks:
		fig = plt.figure(figsize=(10,9))
		fig.suptitle('%s method'%method,fontsize=21)
		plt.loglog(10**freqs, 10**Periodogram, 'k.')
		for frec_indx in result:
			plt.axvline(x=10**freqs[frec_indx-1],color='b')
		plt.xlabel(r'Frequency (h$^{-1}$)')
		plt.ylabel(r'Periodogram (Jy$^2$)')
		plt.savefig('%s_ruptures_%s.png'%(DATname[:-4],method))
		plt.close()
	del algo
	return result

def Periodogram1Slope_and_Normalise(freqs=[],Periodogram=[],DATname='Periodogram.dat', method = 'BottomUp', lowfreq_indx=1, PlotLinFit = False):
	'''
		Function to compute slope of the periodogram
		Parameters:
			- DATname: name of the text file with the data of the periodogram; used for naming the .dat file with the normalised periodogram, and the figs.
			- method: method used to get the break point in the periodogram
			- lowfreq_indx: index where linear fitting begins (to avoid noise at low freq)
			- PlotLinFit: True to plot linear fit + periodogram
		Returns:
			- slope of the periodogram
			- break index (freq. for which noise dominates the signal)
	'''
	### Call function to get 1 break point in the data
	split_index = getBreakIndx(Periodogram,method,11,1,False)[0]
#	break_index = getBreakIndx(Periodogram,method,11,1,False)[0]
#	if freq0 !=0.:
#		split_index = break_index if break_index > (np.argmin(np.abs(freqs-freq0))+1) else (np.argmin(np.abs(freqs-freq0))+3)
#		print('\n\n break index is %i; indx of freq0 is %i; split indx is %i'%(break_index,np.argmin(np.abs(freqs-freq0))+1,split_index))
#		if split_index > (len(freqs)-1):
#			split_index = len(freqs)-1
#	else:
#		split_index = break_index
	
	x1, y1 = freqs[lowfreq_indx:split_index], Periodogram[lowfreq_indx:split_index]
	
	# Perform the "linear" fit for the first segment
	params_linear, covariance_linear = curve_fit(linear_fit, x1, y1)
	slope, intercept = params_linear
	# Calculate the uncertainties (1 sigma)
	sigma_linear = np.sqrt(np.diag(covariance_linear))
	
	### Compute, save and plot normalised periodogram
	NormFreqs = 10**np.array(freqs[:split_index])
	NormPeriodogram = np.array(10**Periodogram[:split_index]/(10**(linear_fit(freqs[:split_index], *params_linear))))
	Outdata = np.column_stack((NormFreqs, NormPeriodogram))
	np.savetxt(DATname[:-4]+'_normalised_%s.dat'%method, Outdata, delimiter='\t', header='#\t Freq (h^-1) \t Periodogram', comments='')
	fig = pl.figure(figsize=(10,7))
	sub1 = fig.add_subplot(111)
	fig.subplots_adjust(wspace=0.01,hspace=0.01,right=0.98,left=0.125)
	fig.suptitle(DATname[:-4]+' normalised %s'%method,fontsize=18)
	sub1.loglog(NormFreqs,NormPeriodogram,'.k')
	sub1.set_xlabel(r'Frequency (h$^{-1}$)')
	sub1.set_ylabel('Norm. Periodogram')
	pl.savefig(DATname[:-4]+'_normalised_%s.png'%method)
	pl.close()
	
	# Plot the data and fits
	if PlotLinFit:
		fig = plt.figure(figsize=(10,9))
		fig.suptitle('%s %s method'%(DATname[:-4],method),fontsize=21)
		fig.subplots_adjust(wspace=0.01,hspace=0.01,right=0.99,left=0.15,top=0.95,bottom=0.2)
		plt.loglog(10**freqs, 10**Periodogram, 'k.', label='Data')
		# plot slope +- 5 sigma
		x1_extend, y1_extend = freqs[lowfreq_indx-1:split_index+1], Periodogram[lowfreq_indx-1:split_index+1]
		plt.loglog(10**x1_extend, 10**linear_fit(x1_extend, *params_linear), label=f'Slope: {slope:.3f} ± {sigma_linear[0]:.3f}, log freq0: {intercept:.3f} ± {sigma_linear[1]:.3f}')
		plt.fill_between(10**x1_extend, 10**linear_fit(x1_extend, slope - sigma_linear[0], intercept - sigma_linear[1]), 10**linear_fit(x1_extend, slope + sigma_linear[0], intercept + sigma_linear[1]), alpha=0.2, color='blue')
		plt.axvline(x=10**freqs[split_index-1],color='b',ls='--')
		plt.xlabel(r'Frequency (h$^{-1}$)')
		plt.ylabel(r'Periodogram (Jy$^2$)')
		fig.legend(loc='lower center', bbox_to_anchor=(0.5, 0))
		plt.savefig('%s_%s.png'%(DATname[:-4],method))
		plt.close()
		del x1_extend,y1_extend
	del split_index,x1,y1,params_linear,covariance_linear,intercept,sigma_linear,NormFreqs,NormPeriodogram,Outdata
	return slope



# Periodic Signal Simulator:
def PeriodicSignal(times,X0,freq0):
	return np.array(X0*np.sin(freq0*times))

def Periodicity_and_NormPeriodogram(DATname='signal.dat', X0_FACT=[3.],PERIOD_LEN=[0.01], method='Kernrbf', select_periodogram='detrended'):
	'''
		Function to simulate a Light Curve
		Parameters:
		- DATname: name of the text file with the data of the simulated signal
		- X0_FACT: the amplitude of the periodic signal will be X0_FACT times the standard deviation of the random simulated signal
		- PERIOD_LEN: the period of the periodic signal (freq0=1/period) will be PERIOD_LEN times the total lenght (i.e. the duration) of the simulated signal
		- PSlope: slope of the periodogram computed from the simulated signal (before injecting periodicity)
		- POrig: value of the periodogram at freq=0, computed from the simulated signal (before injecting periodicity)
		- PBreak: index of the freq. for which noise dominates the signal: will discard freqs. from that indx
		- select_periodogram: indicate which periodogram will be computed: "original" or "detrended"
	'''
	os.chdir(PATH_sampling)
	data = np.loadtxt(DATname)
	times = np.array([ti for ti in data[:,0]])
	signal = np.array([si for si in data[:,1]])
	
	TOT_TIME = np.max(times) - np.min(times)
	freqs = np.array([2.*np.pi*n/TOT_TIME for n in np.arange(1,int((len(times)+1)/2.))])
	
	### require copy needed of selected data column for C++ processing
	Time = np.require(np.copy(times),requirements=['C','A'])
	PFreq = np.require(np.copy(freqs),requirements=['C','A'])
	Ndata = len(signal)
	Nfreqs = len(freqs)
	
	for amp_fact in X0_FACT:
		PATH_X0 = os.path.join(PATH_sampling,'PSignal_%.1f'%amp_fact)
		if not os.path.exists(PATH_X0):
			os.makedirs(PATH_X0)
		os.chdir(PATH_X0)
		if select_periodogram == 'detrended':
			PATH_X0_PERIODOGRAM = os.path.join(PATH_X0,'Detrended_NormPeriodogram')
		else:
			PATH_X0_PERIODOGRAM = os.path.join(PATH_X0,'NormPeriodogram')
		if not os.path.exists(PATH_X0_PERIODOGRAM):
			os.makedirs(PATH_X0_PERIODOGRAM)
		X0 = amp_fact*np.std(signal)
		for period in PERIOD_LEN:
			print('\n Add periodic signal: %s sim%i noise=%i beta=%.2f X0=%.1f T=%.2f\n'%(sampling,nsim,noiselvl,beta,amp_fact,period))
			
			freq0 = 2.*np.pi/(TOT_TIME*period)
			periodic_signal = PeriodicSignal(times,X0,freq0)
			combined_signal = signal + periodic_signal
			
			if method == METHODS[0]:
				with open(DATname[:-4]+'_X0%.1f_T%.2f.dat'%(amp_fact,period), "w+") as fpOut:
					for ti, amp in zip(times,combined_signal):
						fpOut.write("{:.7e}   {:.7e}\n".format(ti,amp))
				fig = pl.figure(figsize=(10,7))
				sub1 = fig.add_subplot(111)
				fig.subplots_adjust(wspace=0.01,hspace=0.01,right=0.98,left=0.125)
				fig.suptitle('%s sim%i noise=%i beta=%.2f periodicity X0=%.1f T=%.2f'%(sampling,nsim,noiselvl,beta,amp_fact,period),fontsize=15)
				sub1.plot(times,combined_signal,'.k')
				sub1.set_xlabel('Times (h)')
				sub1.set_ylabel('Amplitude')
				pl.savefig(DATname[:-4]+'_Xo%.1f_T%.2f.png'%(amp_fact,period))
				pl.close()
			
			### require copy needed of selected data column for C++ processing
			LCurve = np.require(np.copy(combined_signal),requirements=['C','A'])
			LCurve_err = np.zeros(len(combined_signal)) ### careful: no error used!
			Periodogram = np.require(np.zeros(Nfreqs,dtype=np.float64),requirements=['C','A'])
			### Get periodogram
			if select_periodogram == 'detrended':
				detrendPeriodogramcpp.detrendedPeriodogram(Time,LCurve,LCurve_err, PFreq,Periodogram, Ndata,Nfreqs, 0, 0)
			elif select_periodogram == 'original':
				Periodogramcpp.ComputePeriodogram(Time,LCurve, PFreq,Periodogram, Ndata,Nfreqs)
			else:
				print("ERROR: bad naming of version of the periodogram to compute")
				sys.exit()
			FitFreqs = np.log10(PFreq/(2.*np.pi))
			FitPeriodogram = np.log10(Periodogram)
			### Fit slope, save and plot the normalised periodogram
			os.chdir(PATH_X0_PERIODOGRAM)
			if select_periodogram == 'detrended':
				pdatnam = DATname[:-4]+'_Xo%.1f_T%.2f_periodogram_detrended_%s.dat'%(amp_fact,period,method)
			else:
				pdatnam = DATname[:-4]+'_Xo%.1f_T%.2f_periodogram_%s.dat'%(amp_fact,period,method)
			slope = Periodogram1Slope_and_Normalise(FitFreqs,FitPeriodogram,pdatnam,method,1,PLOT_FIT)
			os.chdir(PATH0)
			with open("periodic_signal_periodogram_slopes.txt", "a") as file:
				file.write("%s\t%s\t%i\t%s\t%.1f\t%.2f\t%s\t%.3f\t%.3f\t%s\n"%(sampling,beta,nsim,noiselvl,amp_fact,period,method,slope,(beta-slope)/beta,select_periodogram))
			os.chdir(PATH_X0)
	os.chdir(PATH_sampling)
	del data,times,signal,TOT_TIME,freqs,Time,PFreq,Ndata,Nfreqs,X0,amp_fact,period,freq0,periodic_signal,combined_signal,LCurve,LCurve_err,Periodogram,FitFreqs,FitPeriodogram,slope

import time
tiempo_inicio = time.time()



with open("periodogram_slopes.txt", "w") as file:
	file.write("Sampling\tbeta\tNsim\tNoise\tMethod\tSlope\tErr.\tVersion\n")

with open("periodic_signal_periodogram_slopes.txt", "w") as file:
	file.write("Sampling\tbeta\tNsim\tNoise\tX0\tPeriod\tMethod\tSlope\tErr.\tVersion\n")
	
for sampling in SAMPLINGS:
	PATH_sampling = os.path.join(PATH0,sampling)
	for beta in BETAS:
		PATH_PERIOD = os.path.join(PATH_sampling,'PERIODOGRAMS')
		for nsim in range(NSIMULATIONS):
			for noiselvl in NOISE_LEVEL:
				for method in METHODS:
					print('\n\n\n Computing periodogram slopes for sampling%s beta=%s SIM%i noiselvl=%s, using the %s method to get the break point\n'%(sampling,beta,nsim,noiselvl,method))
					fname = 'nsampling%s_beta%.2f_sim%i_noise%i_periodogram.dat'%(sampling,beta,nsim,noiselvl)
					os.chdir(PATH_PERIOD)
					data = np.loadtxt(fname) # Load data
					freqs = np.array([np.log10(freq) for freq in data[:,0]])
					Periodogram = np.array([np.log10(pdm) for pdm in data[:,1]])
					slope = Periodogram1Slope_and_Normalise(freqs,Periodogram,fname,method,1,PLOT_FIT)
					os.chdir(PATH0)
					with open("periodogram_slopes.txt", "a") as file:
						file.write("%s\t%s\t%i\t%s\t%s\t%.3f\t%.3f\t%s\n"%(sampling,beta,nsim,noiselvl,method,slope,(beta-slope)/beta,"original"))
					datname = 'nsampling%s_beta%.2f_sim%i_noise%i.dat'%(sampling,beta,nsim,noiselvl)
					Periodicity_and_NormPeriodogram(datname, X0_FACT,PERIOD_LEN, method, 'original')
		
		PATH_PERIOD = os.path.join(PATH_sampling,'DETRENDED_PERIODOGRAMS')
		for nsim in range(NSIMULATIONS):
			for noiselvl in NOISE_LEVEL:
				for method in METHODS:
					print('\n\n\n Computing detrended periodogram slopes for sampling%s beta=%s SIM%i noiselvl=%s, using the %s method to get the break point\n'%(sampling,beta,nsim,noiselvl,method))
					fname = 'nsampling%s_beta%.2f_sim%i_noise%i_periodogram_detrended.dat'%(sampling,beta,nsim,noiselvl)
					os.chdir(PATH_PERIOD)
					data = np.loadtxt(fname) # Load data
					freqs = np.array([np.log10(freq) for freq in data[:,0]])
					Periodogram = np.array([np.log10(pdm) for pdm in data[:,1]])
					slope = Periodogram1Slope_and_Normalise(freqs,Periodogram,fname,method,1,PLOT_FIT)
					os.chdir(PATH0)
					with open("periodogram_slopes.txt", "a") as file:
						file.write("%s\t%s\t%i\t%s\t%s\t%.3f\t%.3f\t%s\n"%(sampling,beta,nsim,noiselvl,method,slope,(beta-slope)/beta,"detrended"))
					datname = 'nsampling%s_beta%.2f_sim%i_noise%i.dat'%(sampling,beta,nsim,noiselvl)
					Periodicity_and_NormPeriodogram(datname, X0_FACT,PERIOD_LEN, method, 'detrended')
					

tiempo_fin = time.time()
tiempo_transcurrido = tiempo_fin - tiempo_inicio
print("El programa ha tardado {:.4f} segundos en ejecutarse.".format(tiempo_transcurrido))


"""
Description of methods:
- method == "Dynamic":
	Detecting Change Points With Dynamic Programming
		DynP leverages a dynamic programming approach to efficiently order the search over all possible segmentations, which helps optimize the process and provide accurate results.
		It works by systematically examining all possible segmentations of a given signal to find the exact minimum of the sum of costs associated with each segmentation.
		The user must specify the number of change points in advance.
- method == "Binary":
	Detecting Change Points With Binary Segmentation
		Binary segmentation is pretty simple: first, it looks for a single change point in the entire signal.
		Once it finds that point, it splits the signal into two parts and repeats the process for each of those parts.
		This keeps going until no more change points are found or a specified stopping criterion is met.
		It has a low complexity, which means it doesn’t take too much time or computing power to run and is a good option for large datasets.
		One downside is that it can sometimes miss change points or detect false ones, especially when the changes are close together or the signal is noisy.
- method == "BottomUp":
	Detecting Change Points With The Bottom-Up Method
		Bottom-up segmentation is an interesting approach for detecting change points in a signal.
		Instead of starting with a single segment and dividing it, bottom-up segmentation begins with numerous change points and gradually reduces them.
		Here’s how it goes: first, the signal is split into smaller pieces along a regular pattern.
		Then, these pieces are progressively combined based on their similarities.
		If two adjacent segments are alike, they’re joined into one larger segment.
		This process continues until we’re left with the most meaningful change points
- method == "Kernlin" or "Kernrbf" or "Kerncos":
	Detecting Change Points With Kernel Change Detection
		Kernel change point detection works by mapping the signal to a higher-dimensional space, called a Reproducing Kernel Hilbert Space (RKHS).
		This might sound fancy, but it’s just a way to make change point detection more powerful and versatile by analyzing the signal from a different perspective.
		The core idea is to find the change points that minimize a certain cost function, which measures how well the segmented signal fits the data.
		When the number of change points is known, you can solve a specific optimization problem to find the best segmentation.
		If you don’t know the number of change points, you can use a penalty term to balance the goodness of fit against the complexity of the segmentation.
		There are three kernel functions available in ruptures: linear, rbf, and cosine.
"""
