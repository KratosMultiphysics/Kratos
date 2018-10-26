//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Anna Rehr
//

#if !defined(KRATOS_FORMFINDING_UTILIY_H_INCLUDED )
#define  KRATOS_FORMFINDING_UTILIY_H_INCLUDED



// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"


namespace Kratos
{
	class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION)FormfindingIOUtility
	{
		public:
		
		/// Pointer definition
    	KRATOS_CLASS_POINTER_DEFINITION(FormfindingIOUtility);

		FormfindingIOUtility(ModelPart& rModelPart, const Parameters rParameter);

		~FormfindingIOUtility(){}

		void PrintModelPart();

		void ReadPrestressData();

		void PrintPrestressData();

		private:

		ModelPart& mModelPart;


	};  // class Formfinding Utilities.

}  // namespace Kratos.

#endif // KRATOS_FORMFINDING_UTILIY_H_INCLUDED  defined


