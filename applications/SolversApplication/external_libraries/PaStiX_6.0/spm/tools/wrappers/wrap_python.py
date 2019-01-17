#!/usr/bin/env python
"""
Wrapper Python
==============

 @file wrappers/wrap_python.py

 PaStiX generator for the python wrapper

 @copyright 2017-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.0
 @author Mathieu Faverge
 @date 2017-05-04

"""
import os
import re
import argparse
from . import *

indent="    "
iindent=4

# translation_table of types
types_dict = {
    "int":            ("c_int"),
    "int8_t":         ("c_int8"),
    "spm_coeftype_t": ("c_int"),
    "spm_dir_t":      ("c_int"),
    "spm_trans_t":    ("c_int"),
    "spm_uplo_t":     ("c_int"),
    "spm_diag_t":     ("c_int"),
    "spm_side_t":     ("c_int"),
    "spm_driver_t":   ("c_int"),
    "spm_fmttype_t":  ("c_int"),
    "spm_layout_t":   ("c_int"),
    "spm_normtype_t": ("c_int"),
    "spm_rhstype_t":  ("c_int"),
    "spm_mtxtype_t":  ("c_int"),
    "spm_int_t":      ("__spm_int__"),
    "spmatrix_t":     ("pyspm_spmatrix_t"),
    "size_t":         ("c_size_t"),
    "char":           ("c_char"),
    "double":         ("c_double"),
    "float":          ("c_float"),
    "void":           ("c_void"),
    "MPI_Comm":       ("c_int"),
    "FILE":           ("c_void"),
}

def iso_c_interface_type(arg, return_value, args_list, args_size):
    """Generate a declaration for a variable in the interface."""

    if (arg[1] == "*" or arg[1] == "**"):
        is_pointer = True
    else:
        is_pointer = False

    f_type = types_dict[arg[0]]
    if is_pointer:
        if f_type == "c_void":
            f_type = "c_void_p"
        elif f_type == "c_char":
            f_type = "c_char_p"
        elif f_type == "c_int":
            f_type = "c_int_p"
        else:
            f_type = "POINTER("+f_type+")"

    if (not return_value and arg[1] != "**"):
        f_pointer = "value"
    else:
        f_pointer = ""

    f_name = format("\"%s\", " % arg[2] )

    args_size[0] = max(args_size[0], len(f_name))
    args_size[1] = max(args_size[1], len(f_type))
    args_list.append( [ f_name, f_type ] );

def iso_c_wrapper_type(arg, args_list, args_size):
    """Generate a declaration for a variable in the Fortran wrapper."""

    if (arg[1] == "*" or arg[1] == "**"):
        is_pointer = True
    else:
        is_pointer = False

    isfile = (arg[0] == "FILE")
    f_type = types_dict[arg[0]]

    if is_pointer:
        if f_type == "c_void":
            f_type = "c_void_p"
        elif f_type == "c_char":
            f_type = "c_char_p"
        elif f_type == "c_int":
            f_type = "c_int_p"
        else:
            f_type = "POINTER("+f_type+")"

    f_name = arg[2]
    f_call = f_name

    if isfile:
        f_name = ""
        f_call = "None"

    # detect array argument
    if (is_pointer and f_name in arrays_names_2D):
        f_call = f_name + ".ctypes.data_as( " + f_type + " )"
    elif (is_pointer and f_name in arrays_names_1D):
        f_call = f_name + ".ctypes.data_as( " + f_type + " )"

    if arg[1] == "**":
        f_call = "pointer( " + f_call + " )"

    args_size[0] = max(args_size[0], len(f_name))
    args_size[1] = max(args_size[1], len(f_type))
    args_size[2] = max(args_size[2], len(f_call))
    args_list.append( [f_name, f_type, f_call ] )

class wrap_python:

    @staticmethod
    def header( f ):
        filename = os.path.basename( f['filename'] )
        filename = re.sub(r"\.in", "", filename)
        header = '''"""

 @file ''' + filename + '''

 ''' + f['description'] + '''

 @copyright 2017-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.0
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Louis Poirel
 @date 2017-05-04

This file has been automatically generated with gen_wrappers.py

"""
from ctypes import *
import numpy as np
'''
        if f['header'] != "":
            header += "\n" + f['header']
        return header;

    @staticmethod
    def footer( f ):
        filename = os.path.basename( f['filename'] )
        modname = re.sub(r".f90", "", filename, flags=re.IGNORECASE)
        footer = f['footer']
        return footer

    @staticmethod
    def enum( f, enum ):
        """Generate an interface for an enum.
           Translate it into constants."""

        ename  = enum[0]
        params = enum[1]

        # initialize a string with the fortran interface
        py_interface = "class " + ename + ":\n"

        # loop over the arguments of the enum to get max param length
        length=0
        for param in params:
            # Convert IPARM/DPARM to lower case
            if param[0][1:5] == "PARM":
                param[0] = re.sub(r"[ID]PARM_", "", param[0])
                param[0] = param[0].lower()

            # Remove Pastix from everything
            param[0] = re.sub(r"Pastix", "", param[0])
            param[0] = re.sub(r"Spm", "", param[0])

            if ename == "error":
                param[0] = re.sub(r"PASTIX_", "", param[0])
                param[0] = re.sub(r"SPM_", "", param[0])
                param[0] = re.sub(r"ERR_", "", param[0])
            elif ename == "fact_mode" or ename == "factotype" or ename == "solv_mode":
                param[0] = re.sub(r"Fact", "", param[0])
                param[0] = re.sub(r"Solv", "", param[0])
                param[0] = re.sub(r"Mode", "", param[0])
            elif ename == "scheduler":
                param[0] = re.sub(r"Sched", "", param[0])
            elif ename == "threadmode":
                param[0] = re.sub(r"Thread", "", param[0])
            elif ename[0:8] == "compress":
                param[0] = re.sub(r"Compress", "", param[0])
                param[0] = re.sub(r"When", "", param[0])
                param[0] = re.sub(r"Method", "", param[0])
            elif ename == "rhstype":
                param[0] = re.sub(r"Rhs", "", param[0])
            elif ename == "trans":
                param[0] = param[0]
            elif ename == "mtxtype":
                param[1] = re.sub(r"Pastix", "trans.", param[1])
                param[1] = re.sub(r"Spm", "trans.", param[1])
            elif ename == "normtype":
                param[0] = re.sub(r"Norm", "", param[0])
            else:
                param[0] = re.sub(ename, "", param[0], flags=re.IGNORECASE)
            length = max( length, len(param[0]))
        fmt="%-"+ str(length) + "s"

        # loop over the arguments of the enum
        for param in params:
            name  = param[0]
            value = str(param[1])

            py_interface += indent + format(fmt % name) + " = " + value + "\n"

        if ename in f['enums']:
            py_interface += f['enums'][ename]

        return py_interface

    @staticmethod
    def struct(struct):
        """Generate an interface for a struct.
           Translate it into a derived type."""

        # initialize a string with the fortran interface
        py_interface = ""

        s = 0
        name = struct[0][2]
        name = re.sub(r"pastix_", "", name)
        py_interface +=  "class pyspm_" + name + "(Structure):\n"

        s = iindent
        py_interface +=  s*" " + "_fields_ = ["
        headline = (s+len("_fields_ = ["))*" "

        slist = []
        ssize = [ 0, 0 ]

        # loop over the arguments of the enum
        for j in range(1,len(struct)):
            iso_c_interface_type(struct[j], True, slist, ssize)

        s += iindent
        fmt = "(%-"+ str(ssize[0]) + "s %-"+ str(ssize[1]) +"s)"

        for j in range(0,len(slist)):
            if (j > 0):
                py_interface += ",\n" + headline

            py_interface += format( fmt % (slist[j][0], slist[j][1]) )

        py_interface += " ]\n"

        return py_interface

    @staticmethod
    def function(function):
        """Generate an interface for a function."""

        return_type    = function[0][0]
        return_pointer = function[0][1]

        # is it a function or a subroutine
        if (return_type == "void"):
            is_function = False
        else:
            is_function = True

        c_symbol = function[0][2]
        if "pastix" in c_symbol:
            libname = "libpastix"
            prefix  = "pypastix_"
        elif "spm" in c_symbol:
            libname = "libspm"
            prefix  = "pyspm_"
        else:
            print("ERROR: function name without pastix nor spm")
            return

        # loop over the arguments to compose the different call lines
        func_line = "def " + prefix + c_symbol + "("
        func_line_length = len(func_line)
        func_line_indent = len(func_line)

        args_line = iindent*" " + libname + "." + c_symbol + ".argtypes = ["
        args_line_length = len(args_line)
        args_line_indent = len(args_line)

        if is_function:
            call_line = iindent*" " + "return " + libname + "." + c_symbol + "("

            py_type = types_dict[return_type]
            if return_pointer == "*":
                if py_type == "c_void":
                    py_type = "c_void_p"
                else:
                    py_type = "POINTER("+py_type+")"
            elif  return_pointer == "**":
                py_type = "c_void_p"

            retv_line = iindent*" " + libname + "." + c_symbol + ".restype = " + py_type
        else:
            call_line = iindent*" " + libname + "." + c_symbol + "("
            retv_line = ""
        call_line_length = len(call_line)
        call_line_indent = len(call_line)

        args_list = []
        args_size = [ 0, 0, 0 ]

        # loop over the arguments to compose the first line
        for j in range(1,len(function)):
            iso_c_wrapper_type(function[j], args_list, args_size)

        for j in range(0, len(args_list)):
            # pointers
            arg_name = args_list[j][0]
            arg_type = args_list[j][1]
            arg_call = args_list[j][2]

            isfile = (len(arg_name) == 0)

            if j > 0:
                if not isfile:
                    func_line += ","
                    func_line_length += 2
                args_line += ","
                args_line_length += 2
                call_line += ","
                call_line_length += 2

            # func_line
            l = len(arg_name)
            if not isfile:
                if ((func_line_length + l) > 78):
                    func_line_length = func_line_indent
                    func_line += "\n" + func_line_indent*" "
                func_line += " " + arg_name
                func_line_length += l

            # args_line
            l = len(arg_type)
            if ((args_line_length + l) > 78):
                args_line_length = args_line_indent
                args_line += "\n" + args_line_indent*" "
            args_line += " " + arg_type
            args_line_length += l

            # call_line
            l = len(arg_call)
            if ((call_line_length + l) > 78):
                call_line_length = call_line_indent
                call_line += "\n" + call_line_indent*" "
            call_line += " " + arg_call
            call_line_length += l

        py_interface  = func_line + " ):\n"
        py_interface += args_line + " ]\n"
        if len(retv_line) > 0:
            py_interface += retv_line + "\n"
        py_interface += call_line + " )\n"

        # # add common header
        # py_interface += indent + 2*indent + "use iso_c_binding\n"
        # # import derived types
        # for derived_type in used_derived_types:
        #     py_interface += indent + 2*indent + "import " + derived_type +"\n"
        # py_interface += indent + 2*indent + "implicit none\n"


        # # add the return value of the function
        # if (is_function):
        #     py_interface +=  indent + 2*indent + iso_c_interface_type(function[0], True) + "_c"
        #     py_interface += "\n"

        # # loop over the arguments to describe them
        # for j in range(1,len(function)):
        #     py_interface += indent + 2*indent + iso_c_interface_type(function[j], False)
        #     py_interface += "\n"

        # if (is_function):
        #     py_interface += indent + indent + "end function\n"
        # else:
        #     py_interface += indent + indent + "end subroutine\n"

        # py_interface += indent + "end interface\n"

        return py_interface

    @staticmethod
    def wrapper(function):
        """Generate a wrapper for a function.
           void functions in C will be called as subroutines,
           functions in C will be turned to subroutines by appending
           the return value as the last argument."""

        # is it a function or a subroutine
        if (function[0][0] == "void"):
            is_function = False
        else:
            is_function = True

        c_symbol = function[0][2]
        f_symbol = c_symbol + "_c"

        if (is_function):
            initial_indent_signature = len(indent + "subroutine " + c_symbol + "(") * " "
            initial_indent_call      = len(indent + indent + "info = " + f_symbol + "(") * " "
        else:
            initial_indent_signature = len(indent + "subroutine " + c_symbol + "(") * " "
            initial_indent_call      = len(indent + indent + "call " + f_symbol + "(") * " "

        # loop over the arguments to compose the first line and call line
        signature_line = ""
        call_line = ""
        double_pointers = []
        for j in range(1,len(function)):
            if (j != 1):
                signature_line += ", "
                call_line += ", "

            # do not make the argument list too long
            if (j%9 == 0):
                call_line      += "&\n" + initial_indent_call
                signature_line += "&\n" + initial_indent_signature

            # pointers
            arg_type    = function[j][0]
            arg_pointer = function[j][1]
            arg_name    = function[j][2]

            signature_line += arg_name
            if (arg_pointer == "**"):
                aux_name = arg_name + "_aux"
                call_line += aux_name
                double_pointers.append(arg_name)
            elif (arg_pointer == "*"):
                call_line += "c_loc(" + arg_name + ")"
            else:
                call_line += arg_name

        contains_derived_types = False
        for arg in function:
            if (arg[0] in derived_types):
                contains_derived_types = True

        # initialize a string with the fortran interface
        f_wrapper = ""
        f_wrapper += indent + "subroutine "
        f_wrapper += c_symbol + "("

        # add the info argument at the end
        f_wrapper += signature_line
        if (is_function):
            if (len(function) > 1):
                f_wrapper += ", "

            return_type = function[0][0]
            return_pointer = function[0][1]
            if (return_type == "int"):
                return_var = "info"
            else:
                return_var = return_variables_dict[return_type]

            f_wrapper += return_var

        f_wrapper += ")\n"

        # add common header
        f_wrapper += indent + indent + "use iso_c_binding\n"
        f_wrapper += indent + indent + "implicit none\n"

        # loop over the arguments to describe them
        for j in range(1,len(function)):
            f_wrapper += indent + indent + iso_c_wrapper_type(function[j]) + "\n"

        # add the return info value of the function
        if (is_function):
            if (function[0][1] == "*"):
                f_target = ", pointer"
            else:
                f_target = ""

            f_wrapper += indent + indent + types_dict[return_type] + ", intent(out)" + f_target + " :: " + return_var + "\n"

        f_wrapper += "\n"

        # loop over potential double pointers and generate auxiliary variables for them
        for double_pointer in double_pointers:
            aux_name = double_pointer + "_aux"
            f_wrapper += indent + indent + "type(c_ptr) :: " + aux_name + "\n"
            f_wrapper += "\n"

        if (is_function):
            f_return = return_var
            f_return += " = "
        else:
            f_return = "call "

        # generate the call to the C function
        if (is_function and return_pointer == "*"):
            f_wrapper += indent + indent + "call c_f_pointer(" + f_symbol + "(" + call_line + "), " + return_var + ")\n"
        else:
            f_wrapper += indent + indent + f_return + f_symbol + "(" + call_line + ")\n"

        # loop over potential double pointers and translate them to Fortran pointers
        for double_pointer in double_pointers:
            aux_name = double_pointer + "_aux"
            f_wrapper += indent + indent + "call c_f_pointer(" + aux_name + ", " + double_pointer + ")\n"

        f_wrapper += indent + "end subroutine\n"

        return f_wrapper
