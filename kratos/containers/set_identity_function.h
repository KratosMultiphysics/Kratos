//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                    
//


#if !defined(KRATOS_SET_IDENTITY_FUNCTION_H_INCLUDED )
#define  KRATOS_SET_IDENTITY_FUNCTION_H_INCLUDED



// System includes
#include <functional>



namespace Kratos
{
	///@}
	///@name Kratos Classes
	///@{

	/// Identity function is for having the object also as key in sets
	/**
	*/
	template<class TDataType> class SetIdentityFunction
	{
	public:
		TDataType& operator()(TDataType& data)
		{
			return data;
		}
		const TDataType& operator()(const TDataType& data) const
		{
			return data;
		}
	};

	///@}


}  // namespace Kratos.

#endif // KRATOS_SET_IDENTITY_FUNCTION_H_INCLUDED  defined 
