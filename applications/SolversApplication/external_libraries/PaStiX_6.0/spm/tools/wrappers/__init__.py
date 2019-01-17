"""
Wrappers
========

 @file wrappers/__init__.py

 PaStiX wrapper generators module intialization

 @copyright 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.0
 @author Mathieu Faverge
 @date 2017-05-04

"""

# translation_table with names of auxiliary variables
return_variables_dict = {
    "double":            ("value"),
    "float":             ("value"),
    "pastix_int_t":      ("value"),
    "pastix_order_t":    ("order"),
    "spm_int_t"   :      ("value"),
    "spmatrix_t":        ("spmo"),
}

# global list used to determine derived types
derived_types = [ 'spmatrix_t', 'spm_int_t', 'pastix_int_t', 'pastix_data_t', 'pastix_order_t']

# name arrays which will be translated to assumed-size arrays, e.g. pA(*)
arrays_names_2D = ["pA", "pB", "pC", "pAB", "pQ", "pX", "pAs"]
arrays_names_1D = ["colptr", "rowptr", "loc2glob", "dofs", "row",
                   "iparm", "dparm", "bindtab", "perm", "invp", "schur_list",
                   "rang", "tree" ]

__all__ = [ 'return_variables_dict', 'derived_types', 'arrays_names_1D', 'arrays_names_2D' ]

from .wrap_python  import *
from .wrap_fortran import *
