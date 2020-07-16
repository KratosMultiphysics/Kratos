class MotionStringGenerator(object):
	"""Auxiliary class to generate TAU motion strings for a translatory x-motion
	and/or a pitching oscillation.
	"""

	def __init__(self, deltaT, pitchDeg, thetaDeg, thetaRate):
		self.deltaT       = deltaT
		self.pitchDeg = pitchDeg
		self.thetaDeg = thetaDeg
		self.thetaRate = thetaRate
		#print 'deltaT = ', deltaT
		#print 'thetaDeg = ', thetaDeg
		#print 'thetaRate = ', thetaRate


	def GetMotionString(self,step):
		self.time = step*self.deltaT
		self.thetaInstant     = self.thetaDeg[step]
		self.pitchFreq    = self.thetaRate[step]

		p     = 0.
		q     = np.deg2rad(self.pitchFreq)
		r     = 0.
		phi   = 0.
		theta = np.deg2rad(self.pitchDeg) + np.deg2rad(self.thetaInstant)
		psi   = 0.
		u     = 0
		v     = 0.
		w     = 0.
		dx    = 0.
		dy    = 0.
		dz    = 0.
		motionString=" ".join(map(str, [p,q,r,phi,theta,psi,u,v,w,dx,dy,dz]))
		#print 'step = ', step
        	#print 'motionString = ', motionString
        	#print 'q = ', self.pitchFreq
        	#print 'theta = ', self.thetaInstant

		return motionString

	def __call__(self, step):
		return self.GetMotionString(step)