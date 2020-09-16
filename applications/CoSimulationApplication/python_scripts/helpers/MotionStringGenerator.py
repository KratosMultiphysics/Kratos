# -*- coding: utf-8 -*-
#*******************************************************************************
# $Revision: $
# $Date: $
# $Author: $
#******************************************************************************/
#
#******************************************************************************/
#*******************************************************************************
# import modules Python, TAU, and initialize the python classes 
#*******************************************************************************
import time
import shutil
import sys
import array
import string
import scipy
import numpy as np 
from distutils.version import StrictVersion
import scipy
from scipy.io import netcdf
import math as m
import os
import glob 
import re 


class MotionStringGenerator(object):
	"""Auxiliary class to generate TAU motion strings for a translatory x-motion
	and/or a pitching oscillation.
	"""

	def __init__(self, deltaT, thetaDeg, thetaRate):
		self.deltaT       = deltaT
		self.thetaDeg = thetaDeg
		self.thetaRate = thetaRate
		

	def GetMotionString(self,step):
		self.time = step*self.deltaT
		self.thetaInstant     = self.thetaDeg[step]
		self.pitchFreq    = self.thetaRate[step]

		p     = np.deg2rad(self.pitchFreq)
		q     = 0.
		r     = 0.
		phi   = np.deg2rad(self.thetaInstant)
		theta = 0.
		psi   = 0.
		u     = 0
		v     = 0.
		w     = 0.
		dx    = 0.
		dy    = 0.
		dz    = 0.
		motionString=" ".join(map(str, [p,q,r,phi,theta,psi,u,v,w,dx,dy,dz]))

		return motionString

	def __call__(self, step):
		return self.GetMotionString(step)

