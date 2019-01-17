#!/usr/bin/env python
"""
 @file gen_wrappers.py

 Python and Fortran 90 wrapper generator for some of the solverstack
 libraries, inspired from the PLASMA-OMP fortran generator.

 @copyright 2016-2017 University of Tennessee, US, University of
                      Manchester, UK. All rights reserved.
 @copyright 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.0
 @author Pierre Ramet
 @author Mathieu Faverge
 @date 2017-05-04

"""
import os
import re
import argparse
import wrappers

description = '''\
Generates Fortran 90 and Python wrappers from the spm header files.'''

help = '''\
----------------------------------------------------------------------
Example uses:

   $SPM_SRC_DIR/gen_wrappers.py

----------------------------------------------------------------------
'''

# ------------------------------------------------------------
# command line options
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=description,
    epilog=help )
parser.add_argument('--prefix',        action='store', help='Prefix for variables in Makefile', default='./')
parser.add_argument('args', nargs='*', action='store', help='Files to process')
opts = parser.parse_args()

# exclude inline functions from the interface
exclude_list = [ "inline", "spmIntSort", "spmIntMSort" ]

def polish_file(whole_file):
    """Preprocessing and cleaning of the header file.
       Do not change the order of the regular expressions !
       Works with a long string."""

    clean_file = whole_file

    # borrowed from cfwrapper.py
    # Remove C comments:
    clean_file = re.sub(r"(?s)/\*.*?\*/", "", clean_file)
    clean_file = re.sub(r"//.*", "", clean_file)
    # Remove C directives (multilines then monoline):
    clean_file = re.sub(r"(?m)^#(.*[\\][\n])+.*?$", "", clean_file)
    clean_file = re.sub("(?m)^#.*$", "", clean_file)
    clean_file = re.sub("(?m)#.*", "", clean_file)
    # Remove TABs and overnumerous spaces:
    clean_file = clean_file.replace("\t", " ")
    clean_file = re.sub("[ ]{2,}", " ", clean_file)
    # Remove extern C statement:
    clean_file = re.sub("(?m)^(extern).*$", "", clean_file)
    # Remove empty lines:
    clean_file = re.sub(r"(?m)^\n$", "", clean_file)
    # Merge structs
    clean_file = re.sub(r"(?m)$", "", clean_file)

    # Merge string into single line
    clean_file = re.sub(r"\n", "", clean_file)

    # Split the line based on ";" and "}"
    clean_file = re.sub(r";", "\n", clean_file)
    clean_file = re.sub(r"}", "}\n", clean_file)

    return clean_file

def preprocess_list(initial_list):
    """Preprocessing and cleaning of the header file.
       Works with a list of strings.
       Produces a new list in which each function, enum or struct
       corresponds to a single item."""

    # merge braces
    list1 = []
    merged_line = ""
    nopen = 0
    inStruct = False
    for line in initial_list:

        if (line.find("struct") > -1):
            inStruct = True

        if (inStruct):
            split_character = ","
        else:
            split_character = ""

        nopen += line.count("{") - line.count("}")
        merged_line += line + split_character

        if (nopen <= 0):
            list1.append(merged_line)
            merged_line = ""
            isOpen   = False
            inStruct = False
            nopen = 0

    # merge structs
    list2 = []
    merged_line = ""
    for line in list1:

        merged_line += line

        if (line.find("struct") == -1):
            list2.append(merged_line)
            merged_line = ""

    # clean orphan braces
    list3 = []
    for line in list2:

        if (line.strip() != "}"):
            list3.append(line)

    #print '\n\n'.join(list3)

    return list3

def parse_triple(string):
    """Parse string of
       type (*)name
       into triple of [type, pointer, name, const]"""

    if "const" in string:
        const=1
        string = re.sub(r"const", "", string)
    else:
        const=0

    parts = string.split()
    if (len(parts) < 2 or len(parts) > 3):
        print("Error: Cannot detect type for ", string)

    type_part = str.strip(parts[0])

    if (len(parts) == 2):
        name_with_pointer = str.strip(parts[1])
        if (name_with_pointer.find("**") > -1):
            pointer_part = "**"
            name_part = name_with_pointer.replace("**", "")
        elif (name_with_pointer.find("*") > -1):
            pointer_part = "*"
            name_part    = name_with_pointer.replace("*", "")
        else:
            pointer_part = ""
            name_part    = name_with_pointer

    elif (len(parts) == 3):
        if (str.strip(parts[1]) == "**"):
            pointer_part = "**"
            name_part    = str.strip(parts[2])
        elif (str.strip(parts[1]) == "*"):
            pointer_part = "*"
            name_part    = str.strip(parts[2])
        else:
            print("Error: Too many parts for ", string)

    name_part = name_part.strip()

    return [type_part, pointer_part, name_part, const]


def parse_enums(preprocessed_list):
    """Each enum will be parsed into a list of its arguments."""

    enum_list = []
    for proto in preprocessed_list:

        # extract the part of the function from the prototype
        fun_parts = proto.split("{")

        split_fun = fun_parts[0].strip().split()
        if len(split_fun) == 0:
            continue

        if ((split_fun[0] == "enum") or
            (split_fun[0] == "typedef" and split_fun[1] == "enum")):


            if split_fun[0] == "enum":
                enumname = split_fun[1]
            else:
                enumname = split_fun[2]
            enumname = re.sub(r"_e$", "", enumname)
            enumname = re.sub(r"^pastix_", "", enumname)
            enumname = re.sub(r"^spm_", "", enumname)

            args_string = fun_parts[1];
            args_string = re.sub(r"}", "", args_string)
            args_list = args_string.split(",")
            params = [];
            for args in args_list:
                args = args.strip();
                if (args != ""):
                    values = args.split("=")

                    name = values[0].strip()
                    if (len(values) > 1):
                        value = values[1].strip()
                    else:
                        if (len(params) > 0):
                            value = params[len(params)-1][1] + 1
                        else:
                            value = 0

                    params.append([name, value])

            enum_list.append([enumname, params])

    return enum_list


def parse_structs(preprocessed_list):
    """Each struct will be parsed into a list of its arguments."""

    struct_list = []
    for proto in preprocessed_list:

        # extract the part of the function from the prototype
        fun_parts = proto.split("{")

        if (fun_parts[0].find("struct") > -1):
            args_string = fun_parts[1]
            parts = args_string.split("}")
            args_string = parts[0].strip()
            args_string = re.sub(r"volatile", "", args_string)
            if (len(parts) > 1):
                name_string = parts[1]
                name_string = re.sub(r"(?m),", "", name_string)
                name_string = name_string.strip()
            else:
                print("Error: Cannot detect name for ", proto)
                name_string = "name_not_recognized"

            args_list = args_string.split(",")
            params = [];
            params.append(["struct","",name_string])
            for arg in args_list:
                if (not (arg == "" or arg == " ")):
                    params.append(parse_triple(arg))

            struct_list.append(params)
            wrappers.derived_types.append(name_string)

    # reorder the list so that only defined types are exported
    goAgain = True
    while (goAgain):
        goAgain = False
        for istruct in range(0,len(struct_list)-1):
            struct = struct_list[istruct]
            for j in range(1,len(struct)-1):
                type_name = struct_list[istruct][j][0]

                if (type_name in wrappers.derived_types):

                    # try to find the name in the registered types
                    definedEarlier = False
                    for jstruct in range(0,istruct):
                        struct2 = struct_list[jstruct]
                        that_name = struct2[0][2]
                        if (that_name == type_name):
                            definedEarlier = True

                    # if not found, try to find it behind
                    if (not definedEarlier):
                        definedLater = False
                        for jstruct in range(istruct+1,len(struct_list)-1):
                            struct2 = struct_list[jstruct]
                            that_name = struct2[0][2]
                            if (that_name == type_name):
                                index = jstruct
                                definedLater = True

                        # swap the entries
                        if (definedLater):
                            print("Swapping " + struct_list[istruct][0][2] + " and " + struct_list[index][0][2])
                            tmp = struct_list[index]
                            struct_list[index] = struct_list[istruct]
                            struct_list[istruct] = tmp
                            goAgain = True
                        else:
                            print("Error: Cannot find a derived type " + type_name + " in imported structs.")

    return struct_list


def parse_prototypes(preprocessed_list):
    """Each prototype will be parsed into a list of its arguments."""

    function_list = []
    for proto in preprocessed_list:

        if (proto.find("(") == -1):
            continue

        # extract the part of the function from the prototype
        fun_parts = proto.split("(")
        fun_def  = str.strip(fun_parts[0])

        exclude_this_function = False
        for exclude in exclude_list:
            if (fun_def.find(exclude) != -1):
                exclude_this_function = True

        if (exclude_this_function):
            continue

        # clean keywords
        fun_def = fun_def.replace("static", "")

        # extract arguments from the prototype and make a list from them
        if (len(fun_parts) > 1):
            fun_args = fun_parts[1]
        else:
            fun_args = ""

        fun_args = fun_args.split(")")[0]
        fun_args = fun_args.replace(";", "")
        fun_args = re.sub(r"volatile", "", fun_args)
        fun_args = fun_args.replace("\n", "")
        fun_args_list = fun_args.split(",")

        # generate argument list
        argument_list = []
        # function itself on the first position
        argument_list.append(parse_triple(fun_def))
        # append arguments
        for arg in fun_args_list:
            if (not (arg == "" or arg == " ")):
                argument_list.append(parse_triple(arg))

        # add it only if there is no duplicity with previous one
        is_function_already_present = False
        fun_name = argument_list[0][2]
        for fun in function_list:
            if (fun_name == fun[0][2]):
                is_function_already_present = True

        if (not is_function_already_present):
            function_list.append(argument_list)

    return function_list

def write_module(f, wrapper, generator, enum_list, struct_list, function_list):
    """Generate a single Fortran module. Its structure will be:
       enums converted to constants
       structs converted to derived types
       interfaces of all C functions
       Fortran wrappers of the C functions"""

    modulefile = open( f['filename'], "w" )

    f_interface = generator.header( f )
    modulefile.write(f_interface + "\n")

    # enums
    if (len(enum_list) > 0):
        for enum in enum_list:
            f_interface = generator.enum( f, enum )
            modulefile.write(f_interface + "\n")

    # derived types
    if (len(struct_list) > 0):
        for struct in struct_list:
            f_interface = generator.struct(struct)
            modulefile.write(f_interface + "\n")

    # functions
    if (len(function_list) > 0):
        for function in function_list:
            f_interface = generator.function(function)
            modulefile.write(f_interface + "\n")

        if wrapper:
            modulefile.write("contains\n\n")
            modulefile.write("  " + "! Wrappers of the C functions.\n")

            for function in function_list:
                f_wrapper = generator.wrapper(function)
                modulefile.write(f_wrapper + "\n")

    footer=generator.footer( f )
    modulefile.write(footer + "\n")

    modulefile.close()

    return f['filename']


enums_python_coeftype='''
    @staticmethod
    def getptype ( dtype ):
        np_dict = {
            np.dtype('float32')    : coeftype.Float,
            np.dtype('float64')    : coeftype.Double,
            np.dtype('complex64')  : coeftype.Complex32,
            np.dtype('complex128') : coeftype.Complex64,
        }
        if dtype in np_dict:
            return np_dict[dtype]
        else:
            return -1

    @staticmethod
    def getctype ( flttype ):
        class c_float_complex(Structure):
            _fields_ = [("r",c_float),("i", c_float)]
        class c_double_complex(Structure):
            _fields_ = [("r",c_double),("i", c_double)]

        np_dict = {
            coeftype.Float     : c_float,
            coeftype.Double    : c_double,
            coeftype.Complex32 : c_float_complex,
            coeftype.Complex64 : c_double_complex
        }
        if flttype in np_dict:
            return np_dict[flttype]
        else:
            return -1

    @staticmethod
    def getnptype ( flttype ):
        np_dict = {
            coeftype.Float     : np.dtype('float32'),
            coeftype.Double    : np.dtype('float64'),
            coeftype.Complex32 : np.dtype('complex64'),
            coeftype.Complex64 : np.dtype('complex128')
        }
        if flttype in np_dict:
            return np_dict[flttype]
        else:
            return -1
'''

enums_fortran_footer='''  integer, parameter :: spm_int_t = SPM_INT_KIND

contains

  function spm_getintsize()
    integer :: spm_getintsize
    spm_getintsize = kind(SPM_INT_KIND)
    return
  end function spm_getintsize
'''

spm_enums = {
    'filename' : [ "include/spm_const.h" ],
    'python'   : { 'filename'    : "wrappers/python/spm/enum.py.in",
                   'description' : "SPM python wrapper to define enums and datatypes",
                   'header'      : "# Start with __ to prevent broadcast to file importing enum\n__spm_int__ = @SPM_PYTHON_INTEGER@\n",
                   'footer'      : "",
                   'enums'       : { 'coeftype' : enums_python_coeftype,
                                     'mtxtype'  : "    SymPosDef = trans.ConjTrans + 1\n    HerPosDef = trans.ConjTrans + 2\n" }
    },
    'fortran'  : { 'filename'    : "wrappers/fortran90/src/spm_enums.F90",
                   'description' : "SPM fortran 90 wrapper to define enums and datatypes",
                   'header'      : "  implicit none\n",
                   'footer'      : enums_fortran_footer,
                   'enums'       : { 'mtxtype'  : "    enumerator :: SpmSymPosDef = SpmConjTrans + 1\n    enumerator :: SpmHerPosDef = SpmConjTrans + 2\n" }
    },
}

spm = {
    'filename' : [ "include/spm.h" ],
    'python'   : { 'filename'    : "wrappers/python/spm/__spm__.py",
                   'description' : "SPM python wrapper",
                   'header'      : "from . import libspm\nfrom .enum import __spm_int__\n",
                   'footer'      : "",
                   'enums'       : {}
    },
    'fortran'  : { 'filename'    : "wrappers/fortran90/src/spmf.f90",
                   'description' : "SPM Fortran 90 wrapper",
                   'header'      : "  use spm_enums\n  implicit none\n",
                   'footer'      : "",
                   'enums'       : {}
    },
}

def main():

    # common cleaned header files
    preprocessed_list = []

    # source header files
    for f in [ spm_enums, spm ]:
        preprocessed_list = []
        for filename in f['filename']:

            # source a header file
            c_header_file = open(filename, 'r').read()

            # clean the string (remove comments, macros, etc.)
            clean_file = polish_file(c_header_file)

            # convert the string to a list of strings
            initial_list = clean_file.split("\n")

            # process the list so that each enum, struct or function
            # are just one item
            nice_list = preprocess_list(initial_list)

            # compose all files into one big list
            preprocessed_list += nice_list

            # register all enums
            enum_list = parse_enums(preprocessed_list)

            # register all structs
            struct_list = parse_structs(preprocessed_list)

            # register all individual functions and their signatures
            function_list = parse_prototypes(preprocessed_list)

        # export the module
        modulefilename = write_module( f['fortran'], True,
                                       wrappers.wrap_fortran,
                                       enum_list, struct_list, function_list)
        print( "Exported file: " + modulefilename )

        modulefilename = write_module( f['python'], False,
                                       wrappers.wrap_python,
                                       enum_list, struct_list, function_list)
        print( "Exported file: " + modulefilename )

# execute the program
main()
