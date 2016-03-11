#include <Python.h>
#include <assert.h>
#include <inttypes.h>
#include "tree.h"

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

static PyObject *pyAddEntry(pyGTFtree_t *self, PyObject *args) {
    GTFtree *t = self->t;
    char *chrom = NULL, *name = NULL;
    uint32_t start, end, labelIdx;
    uint8_t strand;
    long lstrand, lstart, lend, llabelIdx;

    if(!(PyArg_ParseTuple(args, "skkskk", &chrom, &lstart, &lend, &name, &lstrand, &llabelIdx))) {
        PyErr_SetString(PyExc_RuntimeError, "pyAddEntry received an invalid or missing argument!");
        return NULL;
    }

    //Convert all of the longs
    if(lstart >= (uint32_t) -1 || lend >= (uint32_t) -1 || lend <= lstart) {
        PyErr_SetString(PyExc_RuntimeError, "pyAddEntry received invalid bounds!");
        return NULL;
    }
    start = (uint32_t) lstart;
    end = (uint32_t) lend;
    if(lstrand != 0 && lstrand != 1 && lstrand != 3) {
        PyErr_SetString(PyExc_RuntimeError, "pyAddEntry received an invalid strand!");
        return NULL;
    }
    strand = (uint8_t) lstrand;
    if(llabelIdx >= (uint32_t) -1) {
        PyErr_SetString(PyExc_RuntimeError, "pyAddEntry received an invalid label idx (too large)!");
        return NULL;
    }
    labelIdx = (uint32_t) llabelIdx;

    //Actually add the entry
    if(addGTFentry(t, chrom, start, end, strand, name, labelIdx)) {
        PyErr_SetString(PyExc_RuntimeError, "pyAddEntry received an error while inserting an entry!");
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *pyVine2Tree(pyGTFtree_t *self, PyObject *args) {
    GTFtree *t = self->t;
    sortGTF(t);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *pyPrintGTFtree(pyGTFtree_t *self, PyObject *args) {
    GTFtree *t = self->t;
    printGTFtree(t);

    Py_INCREF(Py_None);
    return Py_None;
}

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_tree(void) {
    PyObject *res;
    errno = 0;

    if(PyType_Ready(&pyGTFtree) < 0) return NULL;
    res = PyModule_Create(&treemodule);
    if(!res) return NULL;

    Py_INCREF(&pyGTFtree);
    PyModule_AddObject(res, "pyGTFtree", (PyObject *) &pyGTFtree);

    return res;
}
#else
//Python2 initialization
PyMODINIT_FUNC inittree(void) {
    errno = 0; //Sometimes libpython2.7.so is missing some links...
    if(PyType_Ready(&pyGTFtree) < 0) return;
    Py_InitModule3("tree", treeMethods, "A module for handling GTF files for deepTools");
}
#endif
