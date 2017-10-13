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
//   Date:                $Date: 2013-11-04 12:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(RVE_UTILITIES_ELEMENT_INFO_H_INCLUDED)
#define RVE_UTILITIES_ELEMENT_INFO_H_INCLUDED

#include <limits>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/define.h"
#include "includes/serializer.h"

namespace Kratos
{

namespace RveUtilities
{


class SolidElementInfo
{

public:

	///@name Type Definitions
    ///@{

	typedef int IndexType;

	///@}
	
public:
	
	///@name Life Cycle
    ///@{
	
	SolidElementInfo()
		: mElementID(0)
		, mGaussPointID(0)
	{
	}
	
	SolidElementInfo(IndexType elementID, IndexType gaussPointID)
		: mElementID(elementID)
		, mGaussPointID(gaussPointID)
	{
	}
	
	SolidElementInfo(const SolidElementInfo& other)
		: mElementID(other.mElementID)
		, mGaussPointID(other.mGaussPointID)
	{
	}
	
	///@}
	
public:

	///@name Public Operators
    ///@{

	SolidElementInfo& operator = (const SolidElementInfo& other)
	{
		if(this != &other)
		{
			mElementID = other.mElementID;
			mGaussPointID = other.mGaussPointID;
		}
		return *this;
	}
	
	///@}
	
public:
	
	///@name Public Operations
    ///@{
	
	inline virtual std::string GetStringExtension()const
	{
		std::stringstream ss;
		ss << "E_" << mElementID << "_GP_" << mGaussPointID;
		return ss.str();
	}
	
	inline virtual std::string GetInfo()const
	{
		std::stringstream ss;
		ss << "SolidElementInfo:" << std::endl;
		ss << "\tElementID: " << mElementID << std::endl;
		ss << "\tGaussPointID: " << mGaussPointID << std::endl;
		return ss.str();
	}
	
	///@}
	
public:
	
	///@name Public Access
    ///@{
	
	inline const IndexType GetElementID()const
	{
		return mElementID;
	}
	
	inline void SetElementID(IndexType elementID)
	{
		mElementID = elementID;
	}
	
	inline const IndexType GetGaussPointID()const
	{
		return mGaussPointID;
	}
	
	inline void SetGaussPointID(IndexType gaussPointID)
	{
		mGaussPointID = gaussPointID;
	}
	
	///@}
	
protected:

	///@name Protected member Variables
    ///@{
	
	IndexType mElementID;
	IndexType mGaussPointID;
	
    ///@}

public:
	
	///@name Serialization
    ///@{
	
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
		rSerializer.save("EleID", mElementID);
		rSerializer.save("GpID", mGaussPointID);
    }

    virtual void load(Serializer& rSerializer)
    {
		rSerializer.load("EleID", mElementID);
		rSerializer.load("GpID", mGaussPointID);
    }
	
    ///@}
	
};



class ShellElementInfo : public SolidElementInfo
{

public:
	
	///@name Life Cycle
    ///@{
	
	ShellElementInfo()
		: SolidElementInfo()
		, mPlyID(0)
		, mPlyIntegrationPointID(0)
	{
	}
	
	ShellElementInfo(IndexType elementID, IndexType gaussPointID)
		: SolidElementInfo(elementID, gaussPointID)
		, mPlyID(0)
		, mPlyIntegrationPointID(0)
	{
	}
	
	ShellElementInfo(IndexType elementID, IndexType gaussPointID,
	                        IndexType laminaID, IndexType laminaIntegrationPointID)
		: SolidElementInfo(elementID, gaussPointID)
		, mPlyID(laminaID)
		, mPlyIntegrationPointID(laminaIntegrationPointID)
	{
	}
	
	ShellElementInfo(const ShellElementInfo& other)
		: SolidElementInfo(other)
		, mPlyID(other.mPlyID)
		, mPlyIntegrationPointID(other.mPlyIntegrationPointID)
	{
	}

	///@}
	
public:

	///@name Public  Operators
    ///@{

	ShellElementInfo& operator = (const ShellElementInfo& other)
	{
		if(this != &other)
		{
			mElementID = other.mElementID;
			mGaussPointID = other.mGaussPointID;
			mPlyID = other.mPlyID;
			mPlyIntegrationPointID = other.mPlyIntegrationPointID;
		}
		return *this;
	}
	
	///@}
	
public:

	///@name Public  Operations
    ///@{

    virtual std::string GetStringExtension()const
	{
		std::stringstream ss;
		ss << "E_" << mElementID << "_IP_" << mGaussPointID 
		   << "_PLY_" << mPlyID << "_PLYIP_" << mPlyIntegrationPointID;
		return ss.str();
	}
	
    virtual std::string GetInfo()const
	{
		std::stringstream ss;
		ss << "ShellElementInfo:" << std::endl;
		ss << "\tElementID: " << mElementID << std::endl;
		ss << "\tGaussPointID: " << mGaussPointID << std::endl;
		ss << "\tPlyID: " << mPlyID << std::endl;
		ss << "\tPlyIntegrationPointID: " << mPlyIntegrationPointID << std::endl;
		return ss.str();
	}
	
	///@}
	
public:

	///@name Public  Access
    ///@{

	inline const IndexType GetPlyID()const
	{
		return mPlyID;
	}
	
	inline void SetPlyID(IndexType laminaID)
	{
		mPlyID = laminaID;
	}
	
	inline const IndexType GetPlyIntegrationPointID()const
	{
		return mPlyIntegrationPointID;
	}
	
	inline void SetPlyIntegrationPointID(IndexType laminaIntegrationPointID)
	{
		mPlyIntegrationPointID = laminaIntegrationPointID;
	}
	
	///@}
	
protected:

	///@name Protected member Variables
    ///@{
	
	IndexType mPlyID;
	IndexType mPlyIntegrationPointID;
	
    ///@}

protected:

	///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SolidElementInfo );
		rSerializer.save("PlyID", mPlyID);
		rSerializer.save("PlyIpID", mPlyID);
    }

    void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SolidElementInfo );
		rSerializer.load("PlyID", mPlyID);
		rSerializer.load("PlyIpID", mPlyID);
    }

    ///@}
	
};

///@name Input/Output funcitons
///@{

inline std::ostream & operator << (std::ostream& rOStream, const SolidElementInfo& rThis)
{
	return rOStream << rThis.GetInfo();
}

inline std::ostream & operator << (std::ostream& rOStream, const ShellElementInfo& rThis)
{
	return rOStream << rThis.GetInfo();
}

///@}

}

} // namespace Kratos

#endif // RVE_UTILITIES_ELEMENT_INFO_H_INCLUDED
