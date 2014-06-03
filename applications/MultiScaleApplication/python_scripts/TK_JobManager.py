from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
import inspect
from collections import OrderedDict

def GetVariablesAndDofs():
	
	variables = []
	dofs = []
	dofsWithReactions = []
	
	# Loop through all loaded modules
	for aModuleKey in sys.modules.keys():
		
		# Get all members that are 'free functions'
		aModule = sys.modules[aModuleKey]
		allMembers = inspect.getmembers(aModule, inspect.isfunction)
		
		# Loop through all members
		for aMember in allMembers:
			
			# The variable 'aMember' is a tuple:
			# (function name, function address)
			# use aMembers[0] to compare the name of the fuction,
			# and use aMember[1]() to invoke the function if it's
			# what we're searching for.
			if(aMember[0] == "Variables"):
				variables += aMember[1]() # invoke method
			elif(aMember[0] == "Dofs"):
				dofs += aMember[1]() # invoke method
			elif(aMember[0] == "DofsWithReactions"):
				dofsWithReactions += aMember[1]() # invoke method
	
	variables = OrderedDict.fromkeys(variables).keys()
	dofs = OrderedDict.fromkeys(dofs).keys()
	dofsWithReactions = OrderedDict.fromkeys(dofsWithReactions).keys()
	
	return (variables, dofs, dofsWithReactions)