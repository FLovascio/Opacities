#define PY_SSIZE_T_CLEAN

#ifdef __linux__ // Should compile under linux -will test-
#include <Python.h>
#include "numpy/arrayobject.h"
#elif __APPLE__ // python is a framework on macOS requires different include
#include <Python/Python.h>
#include <numpy/arrayobject.h>
#else // I do not have a windows machine so will not support win32_64 in the forseable future
#error "Unsupported compiler"
#endif

#include <iostream>
#include <cstdlib>
#include "Conductivity.hxx"
#include "FileIO.hxx"
#include "Opacity.hxx"

static PyObject * effective_medium(PyObject *self, PyObject *args){
}

static PyObject * opacity(PyObject *self, PyObject *args){
}

PyMODINIT_FUNC PyInit_opacity(void)
{
    PyObject *m;
    static void *PyOpacity_API[PyOpacity_API_pointers];
    PyObject *c_api_object;

    m = PyModule_Create(&opacitymodule);
    if (m == NULL)
        return NULL;

    /* Initialize the C API pointer array */
    PyOpacity_API[PyOpacity_System_NUM] = (void *)PyOpacity_System;

    /* Create a Capsule containing the API pointer array's address */
    c_api_object = PyCapsule_New((void *)PyOpacity_API, "opacitymodule._C_API", NULL);

    if (PyModule_AddObject(m, "_C_API", c_api_object) < 0) {
        Py_XDECREF(c_api_object);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}

PyMODINIT_FUNC PyInit_effective_medium(void)
{
    PyObject *m;
    static void *PyEffectiveMedium_API[PyEffectiveMedium_API_pointers];
    PyObject *c_api_object;

    m = PyModule_Create(&effective_mediummodule);
    if (m == NULL)
        return NULL;

    /* Initialize the C API pointer array */
    PyOpacity_API[PyEffectiveMedium_System_NUM] = (void *)PyEffectiveMedium_System;

    /* Create a Capsule containing the API pointer array's address */
    c_api_object = PyCapsule_New((void *)PyEffectiveMedium_API, "effective_medium._C_API", NULL);

    if (PyModule_AddObject(m, "_C_API", c_api_object) < 0) {
        Py_XDECREF(c_api_object);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}