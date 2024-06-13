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
#include <random>

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
        "Generate Simulated Signal."; /* main docstring */
static char SimSignalcpp_docstring[] =
        "Generate Simulated Signal."; /* function docstring */

/* Available functions */
static PyObject *SimSignalcpp(PyObject *self, PyObject *args);

/* Module specification */
static PyMethodDef module_methods[] = {
        {"SimSignalcpp", SimSignalcpp, METH_VARARGS, SimSignalcpp_docstring},
        {NULL, NULL, 0, NULL}   /* terminated by list of NULLs, apparently */
};


/* Initialize the module */
#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef pc_module_def = {
        PyModuleDef_HEAD_INIT,
        "_LCurveSimulator",     /* m_name */
        module_docstring,       /* m_doc */
        -1,                     /* m_size */
        module_methods,         /* m_methods */
        NULL,NULL,NULL,NULL     /* m_reload, m_traverse, m_clear, m_free */
};
PyMODINIT_FUNC PyInit__LCurveSimulator(void) /* init python module */
{
    PyObject *m = PyModule_Create(&pc_module_def);
    import_array();
    return(m);
}
#else
PyMODINIT_FUNC init_XPCal(void)
{
    import_array();
    PyObject *m = Py_InitModule3("_LCurveSimulator", module_methods, module_docstring);
    if (m == NULL)
        return;
}
#endif


// Trick to allow type promotion below
template <typename T>
struct identity_t { typedef T type; };

//////////////////////////////////
// Step 1: Choose Power Law Spectrum
double PowerSpectrum(double freq, double beta) {
	return pow(2 * M_PI * freq, beta);
};
// MAIN FUNCTION:
static PyObject *SimSignalcpp(PyObject *self, PyObject *args)
{
	// All function variables, except function arguments
	double *times,*freqs, *FTreal, *FTimag, *signal;
	double Sfreq;
	//FILE *fout;

	// Create a random number generator
	std::random_device rd;
	std::mt19937 generator(rd());
	std::normal_distribution<double> randomNormal(0.0, 1.0);
	
	// Python Object to return:
	PyObject *ret;
	ret = Py_BuildValue("i",-1);

	// Function arguments:
	// Inputs:
	PyObject *Pytimes, *Pyfreqs; // Python Objects: O
	int Ntimes, Nfreqs; // integer: i
	float beta;
	// Outputs:
	PyObject *PyFTreal, *PyFTimag, *Pysignal; // Python Objects: O
	// Get arguments from python
	if (!PyArg_ParseTuple(args, "OOiifOOO", &Pytimes, &Pyfreqs, &Ntimes, &Nfreqs, &beta, &PyFTreal, &PyFTimag, &Pysignal))
	{
		printf("FAILED LocNormDCF! Wrong arguments!\n"); fflush(stdout);  return ret;
	};
	// Make sure all python integer, float and double variables are called to C++ 
	printf("Call C++ function with: %i times, %i freqs\n", Ntimes, Nfreqs);
	// Convert python objects to C++ vectors
	times = (double*) PyArray_DATA(Pytimes);
	freqs = (double*) PyArray_DATA(Pyfreqs);
	
	// Outputs: link python objects to C++ vectors
	FTreal = (double*) PyArray_DATA(PyFTreal);
	FTimag = (double*) PyArray_DATA(PyFTimag);
	signal = (double*) PyArray_DATA(Pysignal);
	
	for (int i = 0; i < Nfreqs; ++i) {
		// Step 1: choose power law spectra
		Sfreq = PowerSpectrum(freqs[i],beta);
		// Step 2: get Fourier Transfrom of desired data:
		if (Nfreqs % 2 == 0) { // if even number of data points, Fourier component always real for symmetry reasons
			FTreal[i] = randomNormal(generator)*sqrt(0.5*Sfreq);
		} else {
			FTreal[i] = randomNormal(generator)*sqrt(0.5*Sfreq);
			FTimag[i] = randomNormal(generator)*sqrt(0.5*Sfreq);
		};
		
	};
	// Step 3:  obtain a simulated signal by backward Fourier transform from freq to time domain (compute manually)
	for (int i = 0; i < Ntimes; ++i) {
		double ti = times[i];
		double signal_real=0.;//, signal_imag=0.;
		//if (i==0){
		//	fout = fopen("test.dat","w+");
		//};
		for (int j = 0; j < Nfreqs; ++j) {
			signal_real += (1.0/(2.0*M_PI)) * ( FTreal[j] * cos(2.0*M_PI * freqs[j] * ti) - FTimag[j] * sin(2.0 * M_PI * freqs[j] * ti) );
			signal_real += (1.0/(2.0*M_PI)) * ( FTreal[j] * cos(-2.0*M_PI * freqs[j] * ti) + FTimag[j] * sin(-2.0 * M_PI * freqs[j] * ti) );
			//signal_imag += (1.0/(2.0*M_PI)) * ( FTimag[j] * cos(2.0*M_PI * freqs[j] * ti) + FTreal[j] * sin(2.0 * M_PI * freqs[j] * ti) );
			//signal_imag += (1.0/(2.0*M_PI)) * ( FTimag[j] * cos(-2.0*M_PI * freqs[j] * ti) + FTreal[j] * sin(-2.0 * M_PI * freqs[j] * ti) );
			//if (i==0){ // write in a ".dat" file for testing
			//	fprintf(fout,"time %.4e, freq %.4e, signal %.4e + %.4e i\n",ti,freqs[j],signal_real,signal_imag);
			//}
		};
		signal[i] = signal_real;
		//if (i==0){
		//	fclose(fout);
		//};
	};
	
	ret = Py_BuildValue("i",0);
	return ret;
};

