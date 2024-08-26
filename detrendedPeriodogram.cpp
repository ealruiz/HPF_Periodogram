#include <Python.h>

#if PY_MAJOR_VERSION >= 3
#define NPY_NO_DEPRECATED_API 0x0
#endif
#include <numpy/npy_common.h>
#include <numpy/arrayobject.h>

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <new>
#include <ctime>
#include <sys/stat.h>
#include <string.h>
#include <dirent.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <float.h>
#include <vector>
#include <algorithm> 

// cribbed from SWIG machinery
#if PY_MAJOR_VERSION >= 3
#define PyClass_Check(obj) PyObject_IsInstance(obj, (PyObject *)&PyType_Type)
#define PyInt_Check(x) PyLong_Check(x)
#define PyInt_AsLong(x) PyLong_AsLong(x)
#define PyInt_FromLong(x) PyLong_FromLong(x)
#define PyInt_FromSize_t(x) PyLong_FromSize_t(x)
#define PyString_Check(name) PyBytes_Check(name)
#define PyString_FromString(x) PyUnicode_FromString(x)
#define PyString_Format(fmt, args)  PyUnicode_Format(fmt, args)
#define PyString_Size(str) PyBytes_Size(str)
#define PyString_InternFromString(key) PyUnicode_InternFromString(key)
#define Py_TPFLAGS_HAVE_CLASS Py_TPFLAGS_BASETYPE
#define PyString_AS_STRING(x) PyUnicode_AS_STRING(x)
#define _PyLong_FromSsize_t(x) PyLong_FromSsize_t(x)
#endif

// and after some hacking
#if PY_MAJOR_VERSION >= 3
#define PyString_AsString(obj) PyUnicode_AsUTF8(obj)
#endif

/* Docstrings */
static char module_docstring[] =
        "Process detrended periodogram."; /* main docstring */
static char detrendedPeriodogram_docstring[] =
        "Process detrended periodogram."; /* function docstring */

/* Available functions */
static PyObject *detrendedPeriodogram(PyObject *self, PyObject *args);

/* Module specification */
static PyMethodDef module_methods[] = {
        {"detrendedPeriodogram", detrendedPeriodogram, METH_VARARGS, detrendedPeriodogram_docstring},
        {NULL, NULL, 0, NULL}   /* terminated by list of NULLs, apparently */
};


/* Initialize the module */
#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef pc_module_def = {
        PyModuleDef_HEAD_INIT,
        "_detrendedPeriodogram", /* m_name */
        module_docstring,        /* m_doc */
        -1,                      /* m_size */
        module_methods,          /* m_methods */
        NULL,NULL,NULL,NULL      /* m_reload, m_traverse, m_clear, m_free */
};
PyMODINIT_FUNC PyInit__detrendedPeriodogram(void) /* init python module */
{
    PyObject *m = PyModule_Create(&pc_module_def);
    import_array();
    return(m);
}
#else
PyMODINIT_FUNC init_XPCal(void)
{
    import_array();
    PyObject *m = Py_InitModule3("_detrendedPeriodogram", module_methods, module_docstring);
    if (m == NULL)
        return;
}
#endif


// Trick to allow type promotion below
template <typename T>
struct identity_t { typedef T type; };

//////////////////////////////////
// Calculate second derivatives at the given data points for spline interpolation
std::vector<double> spline_coefficients(const std::vector<double>& bin_times, const std::vector<double>& bin_signal)
{
	int Ndat = bin_signal.size();
	std::vector<double> spc(Ndat, 0.0);
	std::vector<double> u(Ndat-1, 0.0);

	// spc[0] and u[0] are 0.! If not, then use (x - bin_times; y - bin_signal):
	// spc[0] = -0.5;
	// u[0]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	for (int i = 1; i < Ndat - 1; ++i)
	{
		double sig = (bin_times[i]-bin_times[i-1]) / (bin_times[i+1]-bin_times[i-1]);
		double p = sig * spc[i-1] + 2.0;
		spc[i] = (sig - 1.0) / p;
		double uterm1 = (bin_signal[i+1]-bin_signal[i])/(bin_times[i+1]-bin_times[i]);
		double uterm2 = (bin_signal[i]-bin_signal[i-1])/(bin_times[i]-bin_times[i-1]);
		if(std::isnan(uterm1))
		{
			uterm1 = 0;
		};
		if(std::isnan(uterm2))
		{
			uterm2 = 0;
		};
		u[i] = (uterm1 - uterm2);
		u[i] = (6.0*u[i]/(bin_times[i + 1]-bin_times[i-1]) - sig*u[i-1]) / p;
	};
	// spc[Ndat-1] is 0.! If not, then use:
	// qn=0.5;
	// un=(3.0/(x[Ndat-1]-x[Ndat-2]))*(ypn-(y[Ndat-1]-y[Ndat-2])/(x[Ndat-1]-x[Ndat-2]));
	spc[Ndat - 1] = 0.0;

	for (int k = Ndat - 2; k >= 0; --k)
	{
		spc[k] = spc[k] * spc[k + 1] + u[k];
	};
	
	return spc;
};

// Interpolate using the spline function to get values at specified times
double spline_eval(double time, const std::vector<double>& bin_times, const std::vector<double>& bin_signal, const std::vector<double>& spc)
{
	int Ndat = bin_signal.size();
	int klo = 0;
	int khi = Ndat - 1;
	
	while (khi - klo > 1)
	{
		int k = (khi + klo) / 2;
		if (bin_times[k] > time)
		{
			khi = k;
		} else {
			klo = k;
		};
	};
	
	double h = bin_times[khi] - bin_times[klo];
	double a = (bin_times[khi] - time) / h;
	double b = (time - bin_times[klo]) / h;
	
	double spline = a*bin_signal[klo] + b*bin_signal[khi] + ( (a*a*a - a) * spc[klo] + (b*b*b - b) * spc[khi] ) * (h*h)/6.0;
	
	return spline;
};

// Function to detrend a signal by binning the data and interpolate using a spline
std::vector<double> detrend_signal(const std::vector<double>& times, const std::vector<double>& signal, const std::vector<double>& error, int Ndata, double DETREND_BIN = 0.0, int logscale = 0, int save_detrend_signal = 0, const std::string& fname = "detrend.dat")
{
	std::vector<double> detrended_signal(Ndata,0);
	//double signal_min = *std::min_element(signal.begin(), signal.end());
	
	// Binning the data, unless the bin is equal to the duration of the experiment
	std::vector<double> bin_times, bin_signal, bin_error;
	if (DETREND_BIN > 0.0)
	{
		// The first point is estimated independently as the average of the points in the interval (x[0], x[0]+bin/2)
		bin_times.push_back(times[0]);
		double temp_signal = signal[0], temp_error = error[0];
		int nsamebin = 1;
		for (int i = 1; i < Ndata; ++i)
		{
			if (times[i] - times[0] < DETREND_BIN / 2.0)
			{
				temp_signal += signal[i];
				temp_error += error[i];
				++nsamebin;
			} else {
				break;
			};
		};
		bin_signal.push_back(temp_signal/nsamebin);
		bin_error.push_back(temp_error/nsamebin);
		
		int k = 0, flag = 0;
		double temp_times = times[0]; // initialize bin arrays (beginning of 2nd bin)
		temp_signal = signal[0];
		temp_error = error[0];
		nsamebin = 1;
		while (flag != 1)
		{ // loop until all data is binned
			for (int i = k + 1; i < Ndata + 1; ++i)
			{
				if (i == Ndata)
				{
					flag = 1; // stop loop - all data binned, except last bin
					break;
				};
				if (times[i] - times[k] < DETREND_BIN) {
					temp_times += times[i]; // increase value in bin, will be averaged
					temp_signal += signal[i];
					temp_error += error[i];
					++nsamebin; // count number of data in the bin
				} else {
					bin_times.push_back(temp_times / nsamebin); // averages bin arrays
					bin_signal.push_back(temp_signal / nsamebin);
					bin_error.push_back(temp_error / nsamebin);
					k = i; // update k to begin at the proper index in the next loop
					temp_times = times[i]; // init. bin arrays (beginning of next bin)
					temp_signal = signal[i];
					temp_error = error[i];
					nsamebin = 1; // reset count of data in the bin
					break;
				};
			};
		};
		
		// The last two points are estimated independently as the average of the points in the interval (x[Ndata-1]-bin/2, x[Ndata-1]) 
		temp_times = times[Ndata-1]; temp_signal = signal[Ndata-1]; temp_error = error[Ndata-1];
		nsamebin = 1;
		for (int i = Ndata - 2; i > 0; --i)
		{
			if (times[Ndata-1] - times[i] < DETREND_BIN / 2.0)
			{
				temp_times += times[i];
				temp_signal += signal[i];
				temp_error += error[i];
				++nsamebin;
			} else {
				break;
			};
		};
		bin_times.push_back(temp_times / nsamebin);
		bin_signal.push_back(temp_signal / nsamebin);
		bin_error.push_back(temp_error / nsamebin);
	
	} else {
		bin_times = times;
		bin_signal = signal;
		bin_error = error;
	};
	
	int Nbins = bin_times.size();
	
	// Refer all times in arrays to 1st time (i.e. begin at 0). Needed for a good interp.
	double time0 = times[0];
	std::vector<double> times_ini0(Ndata,0);
	std::vector<double> bin_times_ini0(Nbins,0);
	for (int jindx = 0; jindx < Ndata; ++jindx)
	{
		times_ini0[jindx] = (times[jindx] - time0);
	};
	for (int jindx = 0; jindx < Nbins; ++jindx)
	{
		bin_times_ini0[jindx] = (bin_times[jindx] - time0);
	};
	
	// Compute spline coefficients and evaluate spline
	std::vector<double> spline_coeff = spline_coefficients(bin_times_ini0, bin_signal);
	std::vector<double> spline_signal(Ndata,0);
	for (int jindx = 0; jindx < Ndata; ++jindx)
	{
		double time = times_ini0[jindx];
		double tspline = spline_eval(time, bin_times_ini0, bin_signal, spline_coeff);
		if (logscale==1)
		{// returns the spline signal to orig. scale
			spline_signal[jindx] = pow(10,tspline);
		} else {
			spline_signal[jindx] = tspline;
		};
	};
	
	// Compute detrended signal
	for (int i = 0; i < Ndata; ++i)
	{
		if (logscale==1)
		{// detrends the signal by division
			detrended_signal[i] = pow(10,signal[i]) / spline_signal[i]; // *signal_min;
		} else {
			detrended_signal[i] = signal[i] - spline_signal[i];
		};
	};

	// Save detrended signal if required
	if (save_detrend_signal == 1)
	{
		std::ofstream outfile(fname);
		outfile << "#\tJDTIME (h)\t Spline signal (Jy)\t Detrended signal (Jy)\n";
		for (int i = 0; i < Ndata; ++i)
		{
			outfile << times[i] << "\t" << spline_signal[i] << "\t" << detrended_signal[i] << "\n";
		};
		outfile.close();
	};

	return detrended_signal;
};

//////////////////////////////////
// MAIN FUNCTION:
static PyObject *detrendedPeriodogram(PyObject *self, PyObject *args)
{
	// Python Object to return:
	PyObject *ret;
	ret = Py_BuildValue("i",-1);

	// Function arguments:
	// Inputs:	
	PyObject *PyTIME, *PySIGNAL, *PyERROR, *PyFREQS; // Python Objects: O
	int NUMDATA, NUMFREQs, logscale, save_detrend_signal; // integer: i
	// Outputs:
	PyObject *PyPERIODOGRAM; // Python Objects: O
	
	// Get arguments from python
	if (!PyArg_ParseTuple(args, "OOOOOiiii", &PyTIME, &PySIGNAL, &PyERROR, &PyFREQS, &PyPERIODOGRAM, &NUMDATA, &NUMFREQs, &logscale, &save_detrend_signal))
	{
		printf("FAILED LocNormDCF! Wrong arguments!\n"); fflush(stdout);  return ret;
	};
	
	// Convert python objects to C++ vectors
	double *TIMES, *SIGNAL, *ERROR, *freqs, *periodogram;
	TIMES = (double*) PyArray_DATA(PyTIME);
	SIGNAL = (double*) PyArray_DATA(PySIGNAL);
	ERROR = (double*) PyArray_DATA(PyERROR);
	freqs = (double*) PyArray_DATA(PyFREQS);
	// Convert arrays to vectors
	std::vector<double> times(TIMES, TIMES + NUMDATA);
	std::vector<double> signal(SIGNAL, SIGNAL + NUMDATA);
	std::vector<double> error(ERROR, ERROR + NUMDATA);
	std::vector<double> FREQS(freqs, freqs + NUMFREQs);
	// Outputs: link python objects to C++ vectors
	periodogram = (double*) PyArray_DATA(PyPERIODOGRAM);
	
	// Compute the periodogram
	for (int indx = 0; indx < NUMFREQs; ++indx)
	{
		if (indx%111 == 0)
		{
			std::cout << "\rDoing freq. " << indx << " of " << FREQS.size();
			std::cout.flush();
		};
		
		double freq = FREQS[indx];
		double detrend_bin = 2.0 * M_PI / freq;

		// Call detrend_signal function
		std::string foutname = "detrended_signal_freq" + std::to_string(indx) + ".dat";
		std::vector<double> detrended_signal = detrend_signal(times, signal, error, NUMDATA, detrend_bin, logscale, save_detrend_signal, foutname);
		
		// substract the average of the signal: must be 0 to compute the periodogram
		std::vector<double> SIGNAL(NUMDATA, 0);
		double average_detrended_signal = 0.0;
		for (int jindx = 0; jindx < NUMDATA; ++jindx)
		{
			average_detrended_signal += detrended_signal[jindx];
		};
		average_detrended_signal /= NUMDATA;
		for (int jindx = 0; jindx < NUMDATA; ++jindx)
		{
			SIGNAL[jindx] = (detrended_signal[jindx] - average_detrended_signal);
		};
		
		double sum_sin_2wt = 0.0, sum_cos_2wt = 0.0;
		for (double tj : times)
		{
			sum_sin_2wt += std::sin(2.0 * freq * tj);
			sum_cos_2wt += std::cos(2.0 * freq * tj);
		};
		double tau = (1.0/(2.0 * freq)) * std::atan(sum_sin_2wt/sum_cos_2wt);
		
		double sum_cos2_wttau = 0.0, sum_sin2_wttau = 0.0;
		for (double tj : times)
		{
			sum_cos2_wttau += ( std::cos(freq*(tj-tau)) * std::cos(freq*(tj-tau)) );
			sum_sin2_wttau += ( std::sin(freq*(tj-tau)) * std::sin(freq*(tj-tau)) );
		};
		
		double sum_Xcos_wttau = 0.0, sum_Xsin_wttau = 0.0;
		for (int jindx = 0; jindx < NUMDATA; ++jindx)
		{
			sum_Xcos_wttau += (SIGNAL[jindx] * std::cos(freq*(times[jindx]-tau)));
			sum_Xsin_wttau += (SIGNAL[jindx] * std::sin(freq*(times[jindx]-tau)));
		};
		
		double PX = 0.5 * ( (sum_Xcos_wttau*sum_Xcos_wttau / sum_cos2_wttau) + (sum_Xsin_wttau*sum_Xsin_wttau / sum_sin2_wttau) );
		
		periodogram[indx] = PX;
	};

	ret = Py_BuildValue("i",0);
	return ret;
};

