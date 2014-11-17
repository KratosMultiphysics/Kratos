from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
 # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#
#
# Utility routines for automatic benchmarking of Kratos examples        #
#
#
#
# Version history:                                                      #
#
# 1.00: Farshid Mossaiby, Initial Version                               #
# 1.10: Farshid Mossaiby, Added the ability of determining build        #
# reference mode                                #
# 1.20: Farshid Mossaiby, Added self cleaning of BenchTemp.txt          #
#
#

import sys
import os
import time
import string
import types
import smtplib
import math
import platform

Header = "KRATOS_BENCHMARK"
#get the platform used and set the pyhton name 

def TypeToString(Var):
    t = type(Var)
    if (t == float):
        return("Float")
    else:
        if (t == int):
            return("Integer")
        else:
            if (t == bytes):
                return("String")
            else:
                if (t == bool):
                    return("Boolean")
                else:
                    return("Unknown")


def TypedCompare(lr, lt):
    
    EmptyMsg = ""
    matching = True
    non_matching = False

     # Check if it has been an empty line
    if (len(lr) == 1 and len(lt) == 1):
        return matching,EmptyMsg

    # Check if it has been an empty line
    if (len(lr) != len(lt)):
        Msg = "lenght of lines does not match!!\n" 
        Msg += "Reference:" + str(lr) + "\n"
        Msg += "Test:     " + str(lt)
        print(Msg)
        return non_matching,Msg


    
    ref_header = lr[0]
    ref_type = lr[1]
    ref_label = lr[2]
    ref_data = lr[3]
    ref_absolute_tol = lr[4]
    ref_relative_tol = lr[5]

    test_header = lt[0]
    test_type = lt[1]
    test_label = lt[2]
    test_data = lt[3]
    
    try:
        test_absolute_tol = lt[4]
    except:
        test_absolute_tol = "None"  

    try:
        test_relative_tol = lt[5]
    except:
        test_relative_tol = "None"
    
    
    if (ref_type == "Integer"):
        ref_data = int(ref_data)
        test_data = int(test_data)
    else:
        ref_data = float(ref_data)
        test_data = float(test_data)
    


    # Invalid header in reference data?
    if (ref_header != Header.strip()):
        Msg = "Invalid header found for the reference data!\n"
        Msg += "Reference: " + str(lr)
        print(Msg)
        return non_matching,Msg

    # Invalid header in test data?
    if (test_header != Header.strip()):
        Msg = "Invalid header found for the test data!\n"
        Msg += "Test:     " + str(lt)
        print(Msg)
        return non_matching,Msg

    # Incompatible type?
    if (ref_type != test_type):
        Msg = "Incompatible types in benchmark data!\n"
        Msg += "Reference:" + str(lr) + "\n"
        Msg += "Test:     " + str(lt)
        print(Msg)
        return non_matching,Msg

    # Differenet labels?
    if (ref_label != test_label):
        Msg = "Different labels!\n"
        Msg += "Reference:" + str(lr) + "\n"
        Msg += "Test:     " + str(lt)
        print(Msg)
        return non_matching,Msg

    # Differenet tolerances?
    #if (ref_absolute_tol != test_absolute_tol or ref_relative_tol != test_relative_tol):
        #Msg = "Different tolerances!\n"
        #Msg += "Reference:" + str(lr) + "\n"
        #Msg += "Test:     " + str(lt)
        #return non_matching,Msg
    
    exact_types =  ["String","Boolean","Unknown","Integer"]
    if ref_type in exact_types: #do exact check
        if(ref_data == test_data):
            return matching,""
        else:
            Msg = "Difference found in reference and test data:\n"
            Msg += "Reference:" + str(lr) + "\n"
            Msg += "Test:     " + str(lt)
            print(Msg)
            return non_matching,Msg
    else:  #here we need to verify for the tolerance
        
            
        #this shall be the "good way" but it looks like it is too strict!        
        #if(ref_absolute_tol != "None"): ##here we check for identity following comparison guidelines
            ##absolute_tol = float(ref_absolute_tol)
            #epsilon = 1e-6
            #float_MIN_NORMAL = 1e-12 #shouldt define it here, but can not find the equivalent of "limits.h"
            #diff = abs(ref_data - test_data)
            #if (ref_data == test_data): # shortcut, handles infinities
                #return True,""
            #elif(ref_data == 0 or test_data == 0 or diff < float_MIN_NORMAL):
                ## ref_data or test_data is zero or both are extremely close to it
                ## relative error is less meaningful here
                #if(diff < (epsilon * float_MIN_NORMAL)):
                    #return True,""
            #else:
                #if(diff / (abs(ref_data) + abs(test_data) ) < epsilon):
                    #return True,""
        #else:
            ##first of all consider a abosulte tolerance anyway.
            #if(ref_data == 0.0):
                #absolute_tol = 1e-12
            #else:
                #absolute_tol = 1e-12*abs(ref_data)
                
            #if( abs(ref_data - test_data) <= absolute_tol):
                #return matching,""
            
            #if( ref_data != 0.0 and ref_relative_tol != "None"):
                #if(  abs(ref_data - test_data)/abs(ref_data) <= abs(float(ref_relative_tol)) ):
                    #return matching,""

        #first of all consider a abosulte tolerance anyway.
        if(ref_data == 0.0):
            absolute_tol = 1e-12
        else:
            absolute_tol = 1e-12*abs(ref_data)
        
        if(ref_absolute_tol != "None"):
            absolute_tol = float(ref_absolute_tol)
            #print("absolute_tol = ",absolute_tol)
            
            
        if( abs(ref_data - test_data) <= absolute_tol):
            return matching,""
        
        if( ref_data != 0.0 and ref_relative_tol != "None"):
            if(  abs(ref_data - test_data)/abs(ref_data) <= abs(float(ref_relative_tol)) ):
                return matching,""


        Msg = "Difference found in reference and test data:\n"
        Msg += " abs(ref_data - test_data) "+str(abs(ref_data - test_data))+"\n"
        Msg += " absolute_tol "+str(absolute_tol)+"\n"
        Msg += " abs(ref_data - test_data) "+str(abs(ref_data - test_data))+"\n"
        Msg += " abs(ref_data ) "+str(abs(ref_data ))+"\n"
        if(ref_data != 0):
            Msg += " abs(ref_data - test_data)/abs(ref_data ) "+str(abs(ref_data - test_data)/abs(ref_data )) +"\n"
        Msg += " relative_tol "+str(ref_relative_tol)+"\n"      
        Msg += "Reference:" + str(lr) + "\n"
        Msg += "Test:     " + str(lt)
        print(Msg)
        return non_matching,Msg        

    raise Exception("shall never arrive here!!")
        
        #ref_data = lr[3]
        #test_data = lt[3]
        #ref_absolute_tol = lr[4]
        #test_absolute_tol = lt[4]
        #ref_relative_tol = lr[5]
        #test_relative_tol = lt[5]
        
    
    
    
    
    ## Check if it has been an empty line
    #if (len(lr) == 1 and len(lt) == 1):
        #return matching,EmptyMsg

    ## Invalid header in reference data?
    #if (lr[0] != Header.strip()):
        #Msg = "Invalid header found for the reference data!\n"
        #Msg += "Reference: " + str(lr)
        #return non_matching,Msg

    ## Invalid header in test data?
    #if (lt[0] != Header.strip()):
        #Msg = "Invalid header found for the test data!\n"
        #Msg += "Test:     " + str(lt)
        #return non_matching,Msg

    ## Incompatible type?
    #if (lr[1] != lt[1]):
        #Msg = "Incompatible types in benchmark data!\n"
        #Msg += "Reference:" + str(lr) + "\n"
        #Msg += "Test:     " + str(lt)
        #return non_matching,Msg

    ## Differenet labels?
    #if (lr[2] != lt[2]):
        #Msg = "Different labels!\n"
        #Msg += "Reference:" + str(lr) + "\n"
        #Msg += "Test:     " + str(lt)
        #return non_matching,Msg

    ## Differenet tolerances?
    #if (lr[4] != lt[4] or lr[5] != lt[5]):
        #Msg = "Different tolerances!\n"
        #Msg += "Reference:" + str(lr) + "\n"
        #Msg += "Test:     " + str(lt)
        #return non_matching,Msg

    ## Compare based on the type
    
    
    
    ##CHAPUZA! not the way to do it!!
    #if(lr[4] == "None"):
        #lr[4] = 1e-12
        #print("lr[4] =",lr[4])

    ## if comparing strings, booleans or other types, or no tolerance is
    ## specified, do a exact check
    #if (lr[1] == "String" or lr[1] == "Boolean" or lr[1] == "Unknown" or (lr[4] == "None" and lr[5] == "None")):
        #if (lr[3] != lt[3]):
            #Msg = "Difference found in reference and test data:\n"
            #Msg += "Reference:" + str(lr) + "\n"
            #Msg += "Test:     " + str(lt)
            #return non_matching,Msg

    #else:
        #if (lr[1] == "Integer"):
            #m = int(lr[3])
            #n = int(lt[3])
        #else:
            #m = float(lr[3])
            #n = float(lt[3])

        #if (lr[4] == "None"):
            #abserr = True
        #else:
            #abserr = abs(n - m) <= float(lr[4])

        #if (lr[5] == "None"):
            #relerr = True
        #else:
            #if(n != 0.0):
                #relerr = (abs(n - m) / n) <= float(lr[5])
            #else:
                #relerr = True

        #if (not (abserr and relerr)):
            #Msg = "Difference found in reference and test data:\n"
            #Msg += "Reference:" + str(lr) + "\n"
            #Msg += "Test:     " + str(lt)
            #return non_matching,Msg

    #return matching,EmptyMsg


def InBenchmarkingMode():
    return("--benchmarking" in sys.argv)


def InBuildReferenceMode():
    return("--build-reference" in sys.argv)


def Output(Var, Label="", AbsErr=None, RelErr=None):
    if (InBenchmarkingMode()):
        #print(Header, "|", TypeToString(Var), "|", Label, "|", Var, "|", AbsErr, "|", RelErr)

        if(type(Var) == float):
            varout = "{0:.12E}".format(Var)
        else:
            varout = str(Var)
        sys.stdout.flush()
        print(Header, "|", TypeToString(Var), "|", Label, "|", varout, "|", AbsErr, "|", RelErr)
        sys.stdout.flush()
        
def BuildReferenceData(ExamplePath, ReferenceBenchmarkFile):
    if platform.system()=="Windows":
        os.system("python " + ExamplePath + 
        " -u --benchmarking --build-reference | find \"" +
        Header + "\" > " +ReferenceBenchmarkFile)
    else:
        print("BuildReferenceData")        
        if sys.version_info >= (3, 0):
            os.system(
                "python3 " +
                ExamplePath +
                " -u --benchmarking --build-reference | grep \"" +
                Header + "\" > " +
                ReferenceBenchmarkFile)
        else:
            os.system(
                "python -3 " +
                ExamplePath +
                " --benchmarking --build-reference | grep \"" +
                Header + "\" > " +
                ReferenceBenchmarkFile)


def RunBenchmark(ExamplePath, ReferenceBenchmarkFile):
    if platform.system()=="Windows":
        os.system("python " + ExamplePath +
                " -u --benchmarking | find \"" +
                Header +
                "\" > BenchTemp.txt")
    else:
        if sys.version_info >= (3, 0):
            os.system(
                "python3 " +
                ExamplePath +
                " -u --benchmarking | grep \"" +
                Header +
                "\" > BenchTemp.txt")
        else:
            os.system(
                "python -3 " +
                ExamplePath +
                " --benchmarking | grep \"" +
                Header +
                "\" > BenchTemp.txt")
    f = open("BenchTemp.txt", "r")
    t = f.readlines()
    f.close()
    os.remove("BenchTemp.txt")

    f = open(ReferenceBenchmarkFile, "r")
    r = f.readlines()
    f.close()

    t = list(map(str.strip, t))
    r = list(map(str.strip, r))

    successful = True
    if (len(t) != len(r)):
        Msg = "Different amount of benchmark data!"
        successful = False
        return successful,Msg

    n = len(t)

    for i in range(n):
        matching,Msg = TypedCompare(
            list(map(str.strip,
                     r[i].split("|"))),
            list(map(str.strip,
                     t[i].split("|"))))
        if (matching != True):
            successful = False
            return successful,Msg

    Msg = ""
    return successful,Msg


def MPIParallelRunBenchmark(ExamplePath, ReferenceBenchmarkFile):
    if sys.version_info >= (3, 0):
        os.system(
            "mpirun -np 2 python3 " +
            ExamplePath +
            " --benchmarking | grep \"" +
            Header +
            "\" > BenchTemp.txt")
    elif sys.version_info >= (2, 6):
        os.system(
            "mpirun -np 2 python -3 " +
            ExamplePath +
            " --benchmarking | grep \"" +
            Header +
            "\" > BenchTemp.txt")

    f = open("BenchTemp.txt", "r")
    t = f.readlines()
    f.close()
    os.remove("BenchTemp.txt")

    f = open(ReferenceBenchmarkFile, "r")
    r = f.readlines()
    f.close()

    t = list(map(str.strip, t))
    r = list(map(str.strip, r))

    if (len(t) != len(r)):
        Msg = "Different amount of benchmark data!"
        print(Msg)
        return False,Msg

    n = len(t)

    for i in range(n):
        matching,Msg = TypedCompare(
            list(map(str.strip,
                     r[i].split("|"))),
            list(map(str.strip,
                     t[i].split("|"))))
        if (matching != True):
            successful = False
            return False,Msg
    
    return True,""


def StartTiming():
    return (time.time())


def StopTiming(t, AbsDiff, RelDiff):
    if (InBenchmarkingMode()):
        print(Header, "| Time | Timing | ", time.time() - t, "|", AbsDiff, "|", RelDiff)


def NotifyViaEmail(Subject, Text, Recipients):
    print("Sending email, please wait...")
    s = smtplib.SMTP("smtps.cimne.upc.es")
    msg = "From: Kratos benchmarking <no-reply-kratos-benchmarking@cimne.upc.es>\nSubject: " + \
        Subject + "\n" + Text
    s.sendmail(
        "Kratos benchmarking <no-reply-kratos-benchmarking@cimne.upc.es>",
        Recipients,
        msg)
