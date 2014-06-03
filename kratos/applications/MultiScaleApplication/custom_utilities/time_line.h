/*
==============================================================================
KratosMultiScaleApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2014-02-21 10:37:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(TIME_LINE_H_INCLUDED)
#define TIME_LINE_H_INCLUDED



namespace Kratos
{


class AdaptiveTimeLine
{
	
public:
	
	AdaptiveTimeLine()
		: mInitialTime(0.0)
		, mDuration(1.0)
		, mIncrement(1.0)
		, mEndTime(1.0)
		, mCurrentIncrement(1.0)
		, mCurrentTime(0.0)
		, mMinIncrement(1.0)
		, mMaxIncrement(1.0)
		, mTolerance(1.0E-10)
		, mLastIterationConverged(true)
		, mFinished(false)
		, mDone(false)
		, mVerbose(false)
	{
	
	}
	
	AdaptiveTimeLine(double rInitialTime, double rDuration, double rIncrement, 
	                 double rMinIncrement, double rMaxIncrement)
		: mInitialTime(rInitialTime)
		, mDuration(rDuration)
		, mIncrement(rIncrement)
		, mEndTime(0.0)
		, mCurrentIncrement(0.0)
		, mCurrentTime(0.0)
		, mMinIncrement(rMinIncrement)
		, mMaxIncrement(rMaxIncrement)
		, mTolerance(0.0)
		, mLastIterationConverged(true)
		, mFinished(false)
		, mDone(false)
		, mVerbose(false)
	{
		if(mDuration < 0.0) mDuration = -mDuration;
		if(mIncrement < 0.0) mIncrement = -mIncrement;
		
		mEndTime = mInitialTime + mDuration;
		
		mCurrentIncrement = mIncrement;
		mCurrentTime = mInitialTime;
		
		if(mMinIncrement > mMaxIncrement) {
			double temp = mMinIncrement;
			mMinIncrement = mMaxIncrement;
			mMaxIncrement = temp;
		}
		
		mTolerance = mDuration * 1.0E-10;
		if(mTolerance >= mMinIncrement) mTolerance = mMinIncrement * 1.0E-2;
	}
	
public:
	
	inline void SetInitialTime(double rNewInitialTime)
	{
		mInitialTime = rNewInitialTime;
		mEndTime = mInitialTime + mDuration;
		mCurrentTime = mInitialTime;
	}
	
	inline void NextTimeStep(bool rLastIterationConverged)
	{
		mLastIterationConverged = rLastIterationConverged;
		if(mLastIterationConverged)
		{
			mCurrentTime += mCurrentIncrement;
			CheckFinishedState();
			mDone = true;
		}
		else
		{
			if(mVerbose)
			{
				std::stringstream ss;
				ss << " Reducing Time Step due to NON CONVERGENCE: " << std::endl;
				ss << "   Previous Increment : " << mCurrentIncrement << std::endl;
				ss << "   Current Increment  : " << mCurrentIncrement * 0.5 << std::endl;
				std::cout << ss.str();
			}
			mCurrentTime -= mCurrentIncrement;
			mCurrentIncrement *= 0.5;
			if(mCurrentIncrement < mMinIncrement)
			{
				if(mVerbose)
				{
					std::stringstream ss;
					ss << " WARNING: The required increment is smaller than the Mininum increment!" << std::endl;
					ss << " Current increment: " << mCurrentIncrement << " < Min.Increment: " << mMinIncrement << std::endl;
					std::cout << ss.str();
				}
				mDone = false;
			}
			else
			{
				mCurrentTime += mCurrentIncrement;
				CheckFinishedState();
				mDone = true;
			}
		}
	}
	
	inline void NextTimeStep()
	{
		NextTimeStep(true);
	}
	
public:
	
	inline bool Done()const { return mDone; }
	
	inline bool Finished()const { return mFinished; }
	
	inline double CurrentTime()const { return mCurrentTime; }
	
	inline double CurrentIncrement()const { return mCurrentIncrement; }

	inline bool Verbose()const { return mVerbose; }
	
	inline bool& Verbose() { return mVerbose; }
	
private:
	
	inline void CheckFinishedState()
	{
		mFinished = false;
		if(mCurrentTime >= (mEndTime - mTolerance))
		{
			mCurrentTime = mEndTime;
			mFinished = true;
		}
	}
	
private:
	
	double mInitialTime;
	double mDuration;
	double mIncrement;
	double mEndTime;
	double mCurrentIncrement;
	double mCurrentTime;
	double mMinIncrement;
	double mMaxIncrement;
	double mTolerance;
	bool   mLastIterationConverged;
	bool   mFinished;
	bool   mDone;
	bool   mVerbose;
	
};


}


#endif // TIME_LINE_H_INCLUDED