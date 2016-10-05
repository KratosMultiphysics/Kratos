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
//   Date:                $Date: 2013-06-06 10:37:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(LOAD_FUNCTION_H_INCLUDED)
#define LOAD_FUNCTION_H_INCLUDED


#include "boost/python/list.hpp"



namespace Kratos
{


template< class TReal >
class LoadFunction
{

public:

	enum ProlongationType
	{
		Zero,
		Constant,
		Linear
	};

	KRATOS_CLASS_POINTER_DEFINITION( LoadFunction );
	
public:

	LoadFunction() {}
	~LoadFunction() {}

public:

	virtual TReal GetMultiplier(const ProcessInfo& rCurrentProcessInfo)const
	{
		return 1.0;
	}

};


template< class TReal >
class PieceWiseLoadFunction : public LoadFunction< TReal >
{

public:

	KRATOS_CLASS_POINTER_DEFINITION( PieceWiseLoadFunction );

	typedef LoadFunction< TReal > MyBase;
	
	typedef typename LoadFunction< TReal >::ProlongationType ProlongationType;

public:

	PieceWiseLoadFunction(const std::vector<TReal> & Xvalues, const std::vector<TReal> & Yvalues)
		: mX(Xvalues), mY(Yvalues)
		, mProlType(MyBase::Constant)
	{
		CheckInput();
	}
	
	PieceWiseLoadFunction(const std::vector<TReal> & Xvalues, const std::vector<TReal> & Yvalues, ProlongationType prolongationType)
		: mX(Xvalues), mY(Yvalues)
		, mProlType(prolongationType)
	{
		CheckInput();
	}
	
	PieceWiseLoadFunction(const boost::python::list & values)
		: mProlType(MyBase::Constant)
	{
		ProcessPythonList(values);
		CheckInput();
	}

	PieceWiseLoadFunction(const boost::python::list & values, ProlongationType prolongationType)
		: mProlType(prolongationType)
	{
		ProcessPythonList(values);
		CheckInput();
	}

	PieceWiseLoadFunction(const boost::python::list & Xvalues, const boost::python::list & Yvalues, ProlongationType prolongationType)
		: mProlType(prolongationType)
	{
		ProcessPythonList(Xvalues, Yvalues);
		CheckInput();
	}

	~PieceWiseLoadFunction() 
	{
	}

public:

	virtual TReal GetMultiplier(const ProcessInfo& rCurrentProcessInfo)const
	{
		double currentTime = rCurrentProcessInfo[TIME];
		if(currentTime < mX[0]) {
			if(mProlType == MyBase::Zero)
				return 0.0;
			else if(mProlType == MyBase::Constant)
				return mY[0];
			else // linear
				return GetLinearFactorLeft(currentTime);
		}
		else if(currentTime > (mX.back() + mTolerance)) {
			if(mProlType == MyBase::Zero)
				return 0.0;
			else if(mProlType == MyBase::Constant)
				return mY.back();
			else // linear
				return GetLinearFactorRight(currentTime);
		}
		else {
			return GetLinearFactor(currentTime);
		}
	}

private:

	inline void ProcessPythonList(const boost::python::list & pyl)
	{
		size_t n = len(pyl);
		if(mX.size() > 0) mX.clear();
		if(mY.size() > 0) mY.clear();
		if(n < 4) return;
		if(n % 2 != 0) n -= 1;
		n /= 2;
        for(size_t i = 0; i < n; i++) {
			int index = i * 2;
			mX.push_back(boost::python::extract<double>(pyl[index]));
			mY.push_back(boost::python::extract<double>(pyl[index + 1]));
		}
	}
	
	inline void ProcessPythonList(const boost::python::list & plx, const boost::python::list & ply)
	{
        int n = len(plx);
		if(mX.size() > 0) mX.clear();
		if(mY.size() > 0) mY.clear();
		if(n < 2) return;
        if(n != len(ply)) return;
        for(int i = 0; i < n; i++) {
			mX.push_back(boost::python::extract<double>(plx[i]));
			mY.push_back(boost::python::extract<double>(ply[i]));
		}
	}

	inline TReal GetLinearFactor(TReal currentTime)const
	{
		for(size_t i = 1; i < mX.size(); i++) {
			if(currentTime <= mX[i] + mTolerance) {
				TReal x0 = mX[i - 1];
				TReal x1 = mX[i];
				TReal y0 = mY[i - 1];
				TReal y1 = mY[i];
				TReal deltaTime = x1 - x0;
				if(deltaTime == 0.0)
					return y1;
				TReal deltaFactor = y1 - y0;
				TReal relativeTime = currentTime - x0;
				TReal timeRatio = relativeTime / deltaTime;
				return y0 + timeRatio * deltaFactor;
			}
		}
		return 0.0;
	}
	
	inline TReal GetLinearFactorLeft(TReal currentTime)const
	{
		TReal x0 = mX[0];
		TReal x1 = mX[1];
		TReal y0 = mY[0];
		TReal y1 = mY[1];
		TReal deltaTime = x1 - x0;
		if(deltaTime == 0.0)
			return 0.0;
		TReal deltaFactor = y1 - y0;
		TReal relativeTime = currentTime - x0;
		TReal timeRatio = relativeTime / deltaTime;
		return y0 + timeRatio * deltaFactor;
	}
	
	inline TReal GetLinearFactorRight(TReal currentTime)const
	{
		size_t i = mX.size();
		TReal x0 = mX[i - 1];
		TReal x1 = mX[i];
		TReal y0 = mY[i - 1];
		TReal y1 = mY[i];
		TReal deltaTime = x1 - x0;
		if(deltaTime == 0.0)
			return 0.0;
		TReal deltaFactor = y1 - y0;
		TReal relativeTime = currentTime - x0;
		TReal timeRatio = relativeTime / deltaTime;
		return y0 + timeRatio * deltaFactor;
	}

	inline void CheckInput()
	{
		KRATOS_TRY
		
		if(mX.size() != mY.size())
			KRATOS_THROW_ERROR(std::invalid_argument, "PieceWiseLoadFunction - X and Y vectors should have the same size", "");
			
		if(mX.size() < 2)
			KRATOS_THROW_ERROR(std::invalid_argument, "PieceWiseLoadFunction - X vector should have at least 2 values", "");
			
		if(mY.size() < 2)
			KRATOS_THROW_ERROR(std::invalid_argument, "PieceWiseLoadFunction - Y vector should have at least 2 values", "");

		size_t n = mX.size();
		TReal minStep = mX[n - 1] - mX[0];
		for(size_t i = 1; i < mX.size(); i++) {
			if(mX[i] < mX[i - 1])
				KRATOS_THROW_ERROR(std::invalid_argument, "PieceWiseLoadFunction - X vector should be strictly monotonically increasing", "");
			TReal iStep = mX[i] - mX[i - 1];
			if(minStep > iStep)
				minStep = iStep;
		}
		mTolerance = minStep * 1.0E-6;

		KRATOS_CATCH("")
	}
	
private:

	TReal mTolerance;
	std::vector< TReal > mX;
	std::vector< TReal > mY;
	ProlongationType mProlType;
};


}


#endif // LOAD_FUNCTION_H_INCLUDED
