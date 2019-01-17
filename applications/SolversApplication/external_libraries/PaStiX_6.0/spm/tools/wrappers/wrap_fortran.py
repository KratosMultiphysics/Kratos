#!/usr/bin/env python
"""
Wrapper Fortran 90
==================

 @file wrappers/wrap_fortran.py

 PaStiX generator for the Fortran 90 wrapper

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

# set indentation in the f90 file
tab = "  "
indent = "   "

itab=2
iindent=3

# translation_table of types
types_dict = {
    "int":            ("integer(kind=c_int)"),
    "int8_t":         ("integer(kind=c_int8_t)"),
    "spm_coeftype_t": ("integer(c_int)"),
    "spm_dir_t":      ("integer(c_int)"),
    "spm_trans_t":    ("integer(c_int)"),
    "spm_uplo_t":     ("integer(c_int)"),
    "spm_diag_t":     ("integer(c_int)"),
    "spm_side_t":     ("integer(c_int)"),
    "spm_driver_t":   ("integer(c_int)"),
    "spm_fmttype_t":  ("integer(c_int)"),
    "spm_layout_t":   ("integer(c_int)"),
    "spm_normtype_t": ("integer(c_int)"),
    "spm_rhstype_t":  ("integer(c_int)"),
    "spm_mtxtype_t":  ("integer(c_int)"),
    "spmatrix_t":     ("type(spmatrix_t)"),
    "spm_int_t":      ("integer(kind=spm_int_t)"),
    "size_t":         ("integer(kind=c_size_t)"),
    "char":           ("character(kind=c_char)"),
    "double":         ("real(kind=c_double)"),
    "float":          ("real(kind=c_float)"),
    "spm_complex64_t":("complex(kind=c_double_complex)"),
    "spm_complex32_t":("complex(kind=c_float_complex)"),
    "void":           ("type(c_ptr)"),
    "MPI_Comm":       ("integer(kind=c_int)"),
    "FILE":           ("type(c_ptr)"),
}

def iso_c_interface_type(arg, return_value, list):
    """Generate a declaration for a variable in the interface."""

    if (arg[1] == "*" or arg[1] == "**"):
        is_pointer = True
    else:
        is_pointer = False

    if (is_pointer):
        f_type = "type(c_ptr)"
    else:
        f_type = types_dict[arg[0]]

    if (not return_value and arg[1] != "**"):
        f_type += ", "
        f_pointer = "value"
    else:
        f_pointer = ""

    f_name = arg[2]

    list.append( [ f_type, f_pointer, f_name ] );
    return len(f_type + f_pointer)

def iso_c_wrapper_type(arg, args_list, args_size):
    """Generate a declaration for a variable in the Fortran wrapper."""

    if (arg[1] == "*" or arg[1] == "**"):
        is_pointer = True
    else:
        is_pointer = False

    if (is_pointer and (arg[0].strip() == "void")):
        f_type = "type(c_ptr), "
    else:
        f_type = types_dict[arg[0]] + ", "

    if is_pointer and not arg[3]:
        f_intent = "intent(inout)"
    else:
        f_intent = "intent(in)"

    if (is_pointer):
        f_intent += ", "
        if (arg[1] == "*"):
           f_target = "target"
        else:
           f_target = "pointer"
    else:
        f_target = ""

    f_name = arg[2]

    # detect array argument
    if   (is_pointer and f_name in arrays_names_2D):
        f_array = "(:,:)"
    elif (is_pointer and f_name in arrays_names_1D):
        f_array = "(:)"
    else:
        f_array = ""

    f_line = f_type + f_intent + f_target + " :: " + f_name + f_array

    args_size[0] = max(args_size[0], len(f_type))
    args_size[1] = max(args_size[1], len(f_intent))
    args_size[2] = max(args_size[2], len(f_target))
    args_list.append( (f_type, f_intent, f_target, f_name, f_array ) )
    return f_line

class wrap_fortran:

    @staticmethod
    def header( f ):
        filename = os.path.basename( f['filename'] )
        modname = re.sub(r".f90", "", filename, flags=re.IGNORECASE)
        header = '''
!
! @file '''+ filename +'''
!
! ''' + f['description'] + '''
!
! @copyright 2017      Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 6.0.0
! @author Mathieu Faverge
! @date 2017-01-01
!
! This file has been automatically generated with gen_wrappers.py
!
module ''' + modname + '''
  use iso_c_binding
'''

        if f['header'] != "":
            header += f['header']

        return header

    @staticmethod
    def footer( f ):
        filename = os.path.basename( f['filename'] )
        modname = re.sub(r".f90", "", filename, flags=re.IGNORECASE)
        footer = f['footer'] + '''
end module ''' + modname
        return footer

    @staticmethod
    def enum( f, enum ):
        """Generate an interface for an enum.
           Translate it into constants."""

        ename  = enum[0]
        params = enum[1]

        # initialize a string with the fortran interface
        f_interface  = tab + "! enum " + ename + "\n"
        f_interface += tab + "enum, bind(C)\n"

        # loop over the arguments of the enum to get max param length
        length=0
        for param in params:
            length= max( length, len(param[0]) )
        fmt="%-"+ str(length) + "s"

        # Increment for index array enums
        inc = 0
        if ename[1:5] == "parm":
            inc=1

        # loop over the arguments of the enum
        for param in params:
            name  = param[0]
            if isinstance(param[1],int):
                if name[1:10] == "PARM_SIZE":
                    value = str(param[1])
                else:
                    value = str(param[1] + inc)
            else:
                value = str(param[1])
            f_interface += tab + "   enumerator :: " + format(fmt % name) + " = " + value + "\n"

        f_interface += tab + "end enum\n"
        return f_interface

    @staticmethod
    def struct(struct):
        """Generate an interface for a struct.
           Translate it into a derived type."""

        # initialize a string with the fortran interface
        f_interface = ""

        s = itab
        name = struct[0][2]
        f_interface += s*" " + "type, bind(c) :: " + name + "\n"

        # loop over the arguments of the struct to get the length
        s += iindent
        slist = []
        length= 0
        for j in range(1,len(struct)):
            length = max( length, iso_c_interface_type(struct[j], True, slist) )
        fmt = s*" " + "%-"+ str(length) + "s :: %s\n"

        # loop over the arguments of the struct
        for j in range(0,len(slist)):
            f_interface += format( fmt % (slist[j][0], slist[j][2]) )

        s -= iindent
        f_interface += s*" " + "end type " + name + "\n"

        return f_interface

    @staticmethod
    def function(function):
        """Generate an interface for a function."""

        # is it a function or a subroutine
        if (function[0][0] == "void"):
            is_function = False
            symbol="subroutine"
        else:
            is_function = True
            symbol="function"

        c_symbol = function[0][2]
        f_symbol = function[0][2] + "_c"

        used_derived_types = set([])
        for arg in function:
            type_name = arg[0]
            if (type_name in derived_types):
                used_derived_types.add(type_name)

        # initialize a string with the fortran interface
        s = itab
        f_interface = s*" " + "interface\n"
        s += iindent

        f_interface += s*" " + symbol + " " + f_symbol + "("
        s += itab

        initial_len = s + len(symbol + " " + f_symbol + "(" )

        # loop over the arguments to compose the first line
        s += iindent
        for j in range(1,len(function)):
            if (j != 1):
                f_interface += ", "
                initial_len += 2

            l = len(function[j][2])
            if ((initial_len + l) > 77):
                f_interface += "&\n" + s*" "
                initial_len = s

            f_interface += function[j][2]
            initial_len += l

        f_interface += ") &\n"
        f_interface += s*" " + "bind(c, name='" + c_symbol +"')\n"
        s -= iindent

        # add common header
        f_interface += s*" " + "use iso_c_binding\n"
        # import derived types
        for derived_type in used_derived_types:
            f_interface += s*" " + "import " + derived_type +"\n"
        f_interface += s*" " + "implicit none\n"

        plist = []
        length = 0
        # add the return value of the function
        if (is_function):
            l = iso_c_interface_type(function[0], True, plist ) + 2
            length = max( length, l );
            plist[0][2] += "_c"

        # loop over the arguments to describe them
        for j in range(1,len(function)):
            length = max( length, iso_c_interface_type(function[j], False, plist ) )

        # loop over the parameters
        for j in range(0,len(plist)):
            fmt = s*" " + "%s%" + str(length-len(plist[j][0])) + "s :: %s\n"
            f_interface += format( fmt % (plist[j][0], plist[j][1], plist[j][2]) )

        s -= itab
        f_interface += s*" " + "end " + symbol + " " + f_symbol + "\n"

        s -= iindent
        f_interface += s*" " + "end interface\n"

        return f_interface

    @staticmethod
    def wrapper(function):
        """Generate a wrapper for a function.
           void functions in C will be called as subroutines,
           functions in C will be turned to subroutines by appending
           the return value as the last argument."""

        return_type = function[0][0]
        return_pointer = function[0][1]
        return_var = ""

        # is it a function or a subroutine
        if (return_type == "void"):
            is_function = False
            f_return = "call "
        else:
            is_function = True
            if (return_type == "int"):
                return_var = "info"
            else:
                return_var = return_variables_dict[return_type]
            f_return = return_var + " = "

        c_symbol = function[0][2]
        f_symbol = c_symbol + "_c"

        # loop over the arguments to compose the first line and call line
        s = itab
        signature_line = s*" " + "subroutine " + c_symbol + "("
        call_line      = ""
        signature_line_length = len(signature_line)
        call_line_length      = s + itab + len(f_return + f_symbol + "(" )
        double_pointers = []
        args_list = []
        args_size = [ 0, 0, 0 ]
        length = 0
        for j in range(1,len(function)):
            # pointers
            arg_type    = function[j][0]
            arg_pointer = function[j][1]
            arg_name    = function[j][2]

            isfile = (arg_type == "FILE")

            if (j != 1):
                if not isfile:
                    signature_line += ", "
                    signature_line_length += 2
                call_line += ", "
                call_line_length += 2

            # Signature line (do not make the argument list too long)
            if not isfile:
                line = iso_c_wrapper_type(function[j], args_list, args_size)

                l = len(arg_name)
                if ((signature_line_length + l) > 77):
                    signature_line_length = s+itab+iindent
                    signature_line += "&\n" + signature_line_length*" "
                signature_line += arg_name
                signature_line_length += l

            # Call line
            if isfile:
                call_param = "c_null_ptr"
            elif (arg_pointer == "**"):
                aux_name = arg_name + "_aux"
                call_param = aux_name
                double_pointers.append(arg_name)
            elif (arg_pointer == "*") and (not "type(c_ptr), " in args_list[len(args_list)-1]):
                call_param = "c_loc(" + arg_name + ")"
            else:
                call_param = arg_name

            # Do not make the call line too long
            l = len(call_param)
            if ((call_line_length + l) > 77):
                call_line_length = s+itab+iindent+itab
                call_line += "&\n" + call_line_length*" "
            call_line += call_param
            call_line_length += l


        # initialize a string with the fortran interface
        f_wrapper = signature_line
        if (is_function):
            if (len(function) > 1):
                f_wrapper += ", "
                signature_line_length += 2

            return_type = function[0][0]
            return_pointer = function[0][1]

            l = len(return_var)
            if ((signature_line_length + l) > 78):
                signature_line_length = s+itab+iindent
                f_wrapper += "&\n" + signature_line_length*" "

            f_wrapper += return_var

        f_wrapper += ")\n"

        s += itab
        # add common header
        f_wrapper += s*" " + "use iso_c_binding\n"
        f_wrapper += s*" " + "implicit none\n"

        # add the return info value of the function
        if (is_function):
            f_type = types_dict[return_type] + ", "
            f_intent = "intent(out)"
            if (function[0][1] == "*"):
                f_intent += ", "
                f_target = "pointer"
            else:
                f_target = ""

            args_size[0] = max(args_size[0], len(f_type))
            args_size[1] = max(args_size[1], len(f_intent))
            args_size[2] = max(args_size[2], len(f_target))
            args_list.append( (f_type, f_intent, f_target, return_var, "" ) )

        # loop over potential double pointers and generate auxiliary variables for them
        init_line = ""
        for double_pointer in double_pointers:
            aux_name = double_pointer + "_aux"

            args_size[0] = max(args_size[0], len("type(c_ptr)"))
            args_list.append( ("type(c_ptr)", "", "", aux_name, "") )
            init_line += s*" " + aux_name + " = c_loc(" + double_pointer + ")\n"

        # loop over the arguments to describe them
        fmt = s*" " + "%-" + str(args_size[0]) + "s%-" + str(args_size[1]) + "s%-" + str(args_size[2]) + "s :: %s%s\n"
        for j in range(0,len(args_list)):
            f_wrapper += format( fmt % args_list[j] )

        f_wrapper += "\n"
        f_wrapper += init_line
        if init_line != "":
            f_wrapper += "\n"

        # generate the call to the C function
        if (is_function and return_pointer == "*"):
            f_wrapper += s*" " + "call c_f_pointer(" + f_symbol + "(" + call_line + "), " + return_var + ")\n"
        else:
            f_wrapper += s*" " + f_return + f_symbol + "(" + call_line + ")\n"

        # loop over potential double pointers and translate them to Fortran pointers
        for double_pointer in double_pointers:
            aux_name = double_pointer + "_aux"
            f_wrapper += s*" " + "call c_f_pointer(" + aux_name + ", " + double_pointer + ")\n"

        s -= itab
        f_wrapper += s*" " + "end subroutine " + c_symbol + "\n"

        return f_wrapper
