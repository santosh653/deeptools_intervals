#include <Python.h>
#include <structmember.h>
#include "gtf.h"

typedef struct {
    PyObject_HEAD
    GTFtree *t;
} pyGTFtree_t;

/* 
  parseGTF
  parseBED
  findOverlaps

  parseGTF needs to accept the following
    * A file name
    * Two (optional) values (for filtering) in the python interface:
      * "transcript"
      * "exon"

  parseBED needs to intelligently accept 12 columns

  parse* need to handle group labels/multiple files
    * A # line for BED files
    * A deepTools_group field for GTF files

  The python interface to findOverlaps should probably return a tuple.
*/

static void pyGTFDealloc(pyGTFtree_t *self);

static PyMethodDef deeptoolsintervalsMethods[] = {
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
struct deeptoolsintervalsmodule_state {
    PyObject *error;
};

#define GETSTATE(m) ((struct deeptoolsintervalsmodule_state*)PyModule_GetState(m))

static PyModuleDef deeptoolsintervalsmodule = {
    PyModuleDef_HEAD_INIT,
    "deeptools_intervals",
    "A python module creating/accessing GTF-based interval trees with associated meta-data",
    -1,
    deeptoolsintervalsMethods,
    NULL, NULL, NULL, NULL
};
#endif


//Should set tp_dealloc, tp_print, tp_repr, tp_str, tp_members
static PyTypeObject pyGTFtree = {
#if PY_MAJOR_VERSION >= 3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,              /*ob_size*/
#endif
    "deeptoolsintervals.pyGTFtree", /*tp_name*/
    sizeof(pyGTFtree_t),      /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)pyGTFDealloc,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash*/
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    PyObject_GenericGetAttr,   /*tp_getattro*/
    PyObject_GenericSetAttr,   /*tp_setattro*/
    0,                         /*tp_as_buffer*/
#if PY_MAJOR_VERSION >= 3
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
#else
    Py_TPFLAGS_HAVE_CLASS,     /*tp_flags*/
#endif
    "GTF tree",                /*tp_doc*/
    0,                         /*tp_traverse*/
    0,                         /*tp_clear*/
    0,                         /*tp_richcompare*/
    0,                         /*tp_weaklistoffset*/
    0,                         /*tp_iter*/
    0,                         /*tp_iternext*/
    deeptoolsintervalsMethods, /*tp_methods*/
    0,                         /*tp_members*/
    0,                         /*tp_getset*/
    0,                         /*tp_base*/
    0,                         /*tp_dict*/
    0,                         /*tp_descr_get*/
    0,                         /*tp_descr_set*/
    0,                         /*tp_dictoffset*/
    0,                         /*tp_init*/
    0,                         /*tp_alloc*/
    0,                         /*tp_new*/
    0,0,0,0,0,0
};
