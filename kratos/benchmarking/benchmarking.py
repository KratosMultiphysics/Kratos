############################################################################
#                                                                          #
#  Utility routines for automatic benchmarking of Kratos examples          #
#                                                                          #
############################################################################
#                                                                          #
#  Version history:                                                        #
#                                                                          #
#    1.00: Farshid Mossaiby, Initial Version                               #
#    1.10: Farshid Mossaiby, Added the ability of determining build        #
#                            reference mode                                #
#                                                                          #
############################################################################

import sys
import os
import time
import string
import types
import smtplib

Header = "KRATOS_BENCHMARK"

def TypeToString(Var):
	t = type(Var)
	if (t == types.FloatType):
		return("Float")
	else:
		if (t == types.IntType):
			return("Integer")
		else:
			if (t == types.StringType):
				return("String")
			else:
				if (t == types.BooleanType):
					return("Boolean")
				else:
					return("Unknown")

def TypedCompare(lr, lt):
	# Check if it has been an empty line
	if (len(lr) == 1 and len(lt) == 1):
		return True
	
	# Invalid header in reference data?
	if (lr[0] != Header.strip()):
		Msg = "Invalid header found for the reference data!\n"
		Msg += "Reference: " + str(lr)
		return Msg

	# Invalid header in test data?
	if (lt[0] != Header.strip()):
		Msg = "Invalid header found for the test data!\n"
		Msg += "Test:     " + str(lt)
		return Msg

	# Incompatible type?
	if (lr[1] != lt[1]):
		Msg = "Incompatible types in benchmark data!\n"
		Msg += "Reference:" + str(lr) + "\n"
		Msg += "Test:     " + str(lt)
		return Msg

	# Differenet labels?
	if (lr[2] != lt[2]):
		Msg = "Different labels!\n"
		Msg += "Reference:" + str(lr) + "\n"
		Msg += "Test:     " + str(lt)
		return Msg

	# Differenet tolerances?
	if (lr[4] != lt[4] or lr[5] != lt[5]):
		Msg = "Different tolerances!\n"
		Msg += "Reference:" + str(lr) + "\n"
		Msg += "Test:     " + str(lt)
		return Msg

	# Compare based on the type

	# if comparing strings, booleans or other types, or no tolerance is specified, do a exact check
	if (lr[1] == "String" or lr[1] == "Boolean" or lr[1] == "Unknown" or (lr[4] == "None" and lr[5] == "None")):
		if (lr[3] != lt[3]):
			Msg = "Difference found in reference and test data:\n"
			Msg += "Reference:" + str(lr) + "\n"
			Msg += "Test:     " + str(lt)
			return Msg

	else:
		if (lr[1] == "Integer"):
			m = int(lr[3])
			n = int(lt[3])
		else:
			m = float(lr[3])
			n = float(lt[3])

		if (lr[4] == "None"):
			abserr = True
		else:
			abserr = abs(n - m) <= float(lr[4])

		if (lr[5] == "None"):
			relerr = True
		else:
                        if(n != 0.0):
                                relerr = (abs(n - m) / n) <= float(lr[5])
                        else:
                                relerr = True

		if (not (abserr and relerr)):
			Msg = "Difference found in reference and test data:\n"
			Msg += "Reference:" + str(lr) + "\n"
			Msg += "Test:     " + str(lt)
			return Msg

	return True

def InBenchmarkingMode():
	return("--benchmarking" in sys.argv)

def InBuildReferenceMode():
	return("--build-reference" in sys.argv)

def Output(Var, Label = "", AbsErr = None, RelErr = None):
	if (InBenchmarkingMode()):
		print Header, "|", TypeToString(Var), "|", Label, "|", Var, "|", AbsErr, "|", RelErr

def BuildReferenceData(ExamplePath, ReferenceBenchmarkFile):
	os.system("python " + ExamplePath + " --benchmarking --build-reference | grep \"" + Header + "\" > " + ReferenceBenchmarkFile)

def RunBenchmark(ExamplePath, ReferenceBenchmarkFile):
	os.system("python " + ExamplePath + " --benchmarking | grep \"" + Header + "\" > BenchTemp.txt")

	f = open("BenchTemp.txt", "r")
	t = f.readlines()
	f.close()

	f = open(ReferenceBenchmarkFile, "r")
	r = f.readlines()
	f.close()

	t = map(string.strip, t)
	r = map(string.strip, r)

	if (len(t) != len(r)):
		Msg = "Different amount of benchmark data!"
		return Msg

	n = len(t)

	for i in range(n):
		Msg = TypedCompare(map(string.strip, r[i].split("|")), map(string.strip, t[i].split("|")))
		if (Msg != True):
			return Msg
	
	return True

def MPIParallelRunBenchmark(ExamplePath, ReferenceBenchmarkFile):
	os.system("mpirun --np 2 python " + ExamplePath + " --benchmarking | grep \"" + Header + "\" > BenchTemp.txt")

	f = open("BenchTemp.txt", "r")
	t = f.readlines()
	f.close()

	f = open(ReferenceBenchmarkFile, "r")
	r = f.readlines()
	f.close()

	t = map(string.strip, t)
	r = map(string.strip, r)

	if (len(t) != len(r)):
		Msg = "Different amount of benchmark data!"
		return Msg

	n = len(t)

	for i in range(n):
		Msg = TypedCompare(map(string.strip, r[i].split("|")), map(string.strip, t[i].split("|")))
		if (Msg != True):
			return Msg
	
	return True

def StartTiming():
	return (time.time())

def StopTiming(t, AbsDiff, RelDiff):
	if (InBenchmarkingMode()):
		print Header, "| Time | Timing | ", time.time() - t, "|", AbsDiff, "|", RelDiff

def NotifyViaEmail(Subject, Text, Recipients):
	print "Sending email, please wait..."
	s = smtplib.SMTP("smtps.cimne.upc.es")
	msg = "From: Kratos benchmarking <no-reply-kratos-benchmarking@cimne.upc.es>\nSubject: " + Subject + "\n" + Text
	s.sendmail("Kratos benchmarking <no-reply-kratos-benchmarking@cimne.upc.es>", Recipients, msg)

