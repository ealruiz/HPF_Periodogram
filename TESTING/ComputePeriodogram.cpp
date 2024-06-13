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
        "Process periodogram."; /* main docstring */
static char ComputePeriodogram_docstring[] =
        "Process periodogram."; /* function docstring */

/* Available functions */
static PyObject *ComputePeriodogram(PyObject *self, PyObject *args);

/* Module specification */
static PyMethodDef module_methods[] = {
        {"ComputePeriodogram", ComputePeriodogram, METH_VARARGS, ComputePeriodogram_docstring},
        {NULL, NULL, 0, NULL}   /* terminated by list of NULLs, apparently */
};


/* Initialize the module */
#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef pc_module_def = {
        PyModuleDef_HEAD_INIT,
        "_ComputePeriodogram", /* m_name */
        module_docstring,        /* m_doc */
        -1,                      /* m_size */
        module_methods,          /* m_methods */
        NULL,NULL,NULL,NULL      /* m_reload, m_traverse, m_clear, m_free */
};
PyMODINIT_FUNC PyInit__ComputePeriodogram(void) /* init python module */
{
    PyObject *m = PyModule_Create(&pc_module_def);
    import_array();
    return(m);
}
#else
PyMODINIT_FUNC init_XPCal(void)
{
    import_array();
    PyObject *m = Py_InitModule3("_ComputePeriodogram", module_methods, module_docstring);
    if (m == NULL)
        return;
}
#endif


// Trick to allow type promotion below
template <typename T>
struct identity_t { typedef T type; };



//////////////////////////////////
// MAIN FUNCTION:
static PyObject *ComputePeriodogram(PyObject *self, PyObject *args)
{
	// Python Object to return:
	PyObject *ret;
	ret = Py_BuildValue("i",-1);

	// Function arguments:
	// Inputs:	
	PyObject *PyTIME, *PySIGNAL, *PyFREQS; // Python Objects: O
	int NUMDATA, NUMFREQs; // integer: i
	// Outputs:
	PyObject *PyPERIODOGRAM; // Python Objects: O
	
	// Get arguments from python
	if (!PyArg_ParseTuple(args, "OOOOii", &PyTIME, &PySIGNAL, &PyFREQS, &PyPERIODOGRAM, &NUMDATA, &NUMFREQs))
	{
		printf("FAILED LocNormDCF! Wrong arguments!\n"); fflush(stdout);  return ret;
	};
	
	// Convert python objects to C++ vectors
	double *TIMES, *SIGNAL, *freqs, *periodogram;
	TIMES = (double*) PyArray_DATA(PyTIME);
	SIGNAL = (double*) PyArray_DATA(PySIGNAL);
	freqs = (double*) PyArray_DATA(PyFREQS);
	// Convert arrays to vectors
	std::vector<double> times(TIMES, TIMES + NUMDATA);
	std::vector<double> signal(SIGNAL, SIGNAL + NUMDATA);
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

