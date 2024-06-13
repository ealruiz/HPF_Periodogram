import glob, os, gc, re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from fractions import Fraction
import ruptures as rpt
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
DIVIDE_BETA = False # True to generate signals with fixed beta values for the Power Spectrum and divide the signals consequently; if False, will generate signals with an expanded sample of values of beta
if DIVIDE_BETA:
	NSIMULATIONS = 1111 ### number of simulated signals
	BETAS = [-1.5,-2.,-2.5] # fixed amount of betas
else:
	NSIMULATIONS = 333 ### number of simulated signals
	Nbetas = 12 ### will get simulations for Nbetas-1 slopes of the PSD
	BETAS = np.unique(np.round(np.concatenate((np.linspace(-2.5,-2.,int(Nbetas/2)), np.linspace(-2.,-1.5,int(Nbetas/2)))), 2))

NOISE_LEVEL = [3.]

METHODS = ['Dynamic','Binary','BottomUp', 'Kernlin','Kernrbf']

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
#	if method == "Pelt":
#		""" method 1: Detecting Change Points With Pelt
#				Faster method, can use penalty functions, but very simple
#				Caveat for periodogram: too many change points for the periogodram
#		"""
#		algo = rpt.Pelt(model="l2", min_size=MinData) # model to compute distance (l2:rms) and minimum lenght of segment
#		algo.fit(data)
#		result = algo.predict(pen=1) # penality value

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


def Periodogram1Slope(DATname='Periodogram.dat', method = 'BottomUp', lowfreq_indx=1, PlotLinFit = False):
	'''
		Function to compute slope of the periodogram
		Parameters:
			- DATname: name of the text file with the data of the periodogram
			- method: method used to get the break point in the periodogram
			- lowfreq_indx: index where linear fitting begins (to avoid noise at low freq)
			- PlotLinFit: True to plot linear fit + periodogram
		Returns: slope of the periodogram
	'''
	# Load your data
	data = np.loadtxt(DATname)
	
	freqs = np.array([np.log10(freq) for freq in data[:,0]])
	Periodogram = np.array([np.log10(pdm) for pdm in data[:,1]])
	
	### Call function to get 1 break point in the data
	split_index = getBreakIndx(Periodogram,method,11,1,False)[0]
	
	x1, y1 = freqs[lowfreq_indx:split_index], Periodogram[lowfreq_indx:split_index]
	# x1_clean = x1[~np.isnan(y1) & ~np.isinf(y1)]
	# y1_clean = y1[~np.isnan(y1) & ~np.isinf(y1)]
	
	# Perform the "linear" fit for the first segment
	params_linear, covariance_linear = curve_fit(linear_fit, x1, y1)
	slope, intercept = params_linear
	# Calculate the uncertainties (1 sigma)
	sigma_linear = np.sqrt(np.diag(covariance_linear))
	
	# Plot the data and fits
	if PlotLinFit:
		os.chdir(PATH_SLOPE)
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
		os.chdir(PATH_PERIOD)
		del x1_extend,y1_extend
	del data,freqs,Periodogram,split_index,x1,y1,params_linear,covariance_linear,intercept,sigma_linear
	return slope




import time
init_time = time.time()



with open("periodogram_slopes.txt", "w") as file:
	file.write("Sampling\tbeta\tNsim\tNoise\tMethod\tSlope\tErr.\tVersion\n")
	
	for sampling in SAMPLINGS:
		PATH_sampling = os.path.join(PATH0,sampling)
		for beta in BETAS:
			PATH_PERIOD = os.path.join(PATH_sampling,'PERIODOGRAMS')
			os.chdir(PATH_PERIOD)
			PATH_SLOPE = os.path.join(PATH_PERIOD,'SlopeFitting')
			if not os.path.exists(PATH_SLOPE):
				os.makedirs(PATH_SLOPE)
			for nsim in range(NSIMULATIONS):
				for noiselvl in NOISE_LEVEL:
					for method in METHODS:
						print('\nComputing periodogram slopes for sampling%s beta=%s SIM%i noiselvl=%s, using the %s method to get the break point\n'%(sampling,beta,nsim,noiselvl,method))
						fname = 'nsampling%s_beta%.2f_sim%i_noise%i_periodogram.dat'%(sampling,beta,nsim,noiselvl)
						slope = Periodogram1Slope(fname,method,1,False)
						file.write("%s\t%s\t%i\t%s\t%s\t%.3f\t%.3f\t%s\n"%(sampling,beta,nsim,noiselvl,method,slope,(beta-slope)/beta,"original"))
			
			PATH_PERIOD = os.path.join(PATH_sampling,'DETRENDED_PERIODOGRAMS')
			os.chdir(PATH_PERIOD)
			PATH_SLOPE = os.path.join(PATH_PERIOD,'SlopeFitting')
			if not os.path.exists(PATH_SLOPE):
				os.makedirs(PATH_SLOPE)
			for nsim in range(NSIMULATIONS):
				for noiselvl in NOISE_LEVEL:
					for method in METHODS:
						print('\nComputing detrended periodogram slopes for sampling%s beta=%s SIM%i noiselvl=%s, using the %s method to get the break point\n'%(sampling,beta,nsim,noiselvl,method))
						fname = 'nsampling%s_beta%.2f_sim%i_noise%i_periodogram_detrended.dat'%(sampling,beta,nsim,noiselvl)
						slope = Periodogram1Slope(fname,method,1,False)
						file.write("%s\t%s\t%i\t%s\t%s\t%.3f\t%.3f\t%s\n"%(sampling,beta,nsim,noiselvl,method,slope,(beta-slope)/beta,"detrended"))


end_time = time.time()
totat_time_duration = end_time - init_time
print("Total exec. time: {:.4f} seconds.".format(totat_time_duration))


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
- method == "Window":
	Detecting Change Points With Window-Based Method
		Window-based change point detection is a pretty nifty method to find where a signal changes.
		Imagine two windows sliding along the time series, checking out the properties within each window (like the mean).
		What we’re looking for is how similar or different these properties are, and that’s where the discrepancy measure comes in.
		Think of the discrepancy as a way to measure how much difference splitting the signal at a certain point would make.
		If the properties of the windows are similar, the discrepancy will be low. Else, it will be high, and that’s when we suspect a change point occurred.
		One issue is that the choice of window length can have a significant impact on the results.
			Caveat for periodogram: yields worse results
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
