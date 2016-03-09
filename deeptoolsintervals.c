#include <Python.h>
#include <assert.h>
#include <inttypes.h>
#include "deeptoolsintervals.h"

static void pyGTFDealloc(pyGTFtree_t *self) {
    if(self->t) destroyGTFtree(self->t);
    PyObject_DEL(self);
}

#if PY_MAJOR_VERSION >= 3
//Return 1 iff obj is a ready unicode type
int PyString_Check(PyObject *obj) {
    if(PyUnicode_Check(obj)) {
        return PyUnicode_READY(obj)+1;
    }
    return 0;
}

//I don't know what happens if PyBytes_AsString(NULL) is used...
char *PyString_AsString(PyObject *obj) {
    return PyBytes_AsString(PyUnicode_AsASCIIString(obj));
}
#endif

//Will return 1 for long or int types currently
int isNumeric(PyObject *obj) {
#if PY_MAJOR_VERSION < 3
    if(PyInt_Check(obj)) return 1;
#endif
    return PyLong_Check(obj);
}

//On error, throws a runtime error, so use PyErr_Occurred() after this
uint32_t Numeric2Uint(PyObject *obj) {
    long l;
#if PY_MAJOR_VERSION < 3
    if(PyInt_Check(obj)) {
        return (uint32_t) PyInt_AsLong(obj);
    }
#endif
    l = PyLong_AsLong(obj);
    //Check bounds
    if(l > 0xFFFFFFFF) {
        PyErr_SetString(PyExc_RuntimeError, "Length out of bounds for a bigWig file!");
        return (uint32_t) -1;
    }
    return (uint32_t) l;
}

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_deeptoolsintervals(void) {
    PyObject *res;
    errno = 0;

    if(PyType_Ready(&pyGTFtree) < 0) return NULL;
    res = PyModule_Create(&deeptoolsintervalsmodule);
    if(!res) return NULL;

    Py_INCREF(&pyGTFtree);
    PyModule_AddObject(res, "pyGTFtree", (PyObject *) &pyGTFtree);

    return res;
}
#else
//Python2 initialization
PyMODINIT_FUNC initdeeptoolsintervals(void) {
    errno = 0; //Sometimes libpython2.7.so is missing some links...
    if(PyType_Ready(&pyGTFtree) < 0) return;
    Py_InitModule3("deeptoolsintervals", deeptoolsintervalsMethods, "A module for handling GTF files for deepTools");
}
#endif
