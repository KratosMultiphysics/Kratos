//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klauss B Sautter
//                   Philipp Bucher
//


// System includes


// External includes 


// Project includes
#include "static_condensation_utility.h"


namespace Kratos
{
	namespace StaticCondensationUtility
	{
		void CondenseLeftHandSide(){}

        std::vector<MatrixType> CalculateSchurComplements(
			const MatrixType& rLeftHandSideMatrix,
			const std::vector<int> & rDofList)
        {
            
        }

		std::vector<int> CreateRemainingDofList(
			const std::vector<int> & rDofList)
        {

        }

		void FillSchurComplements(
			MatrixType& Submatrix,
			const MatrixType& rLeftHandSideMatrix,
			const std::vector<int>& rVecA,
			const std::vector<int>& rVecB,
			const SizeType& rSizeA,
			const SizeType& rSizeB) //maybe inline
        {

        }

		void ConvertingCondensation(
			Vector& rValues,
			const std::vector<int>& rDofList,
			const MatrixType& rLeftHandSideMatrix)
        {

        }

	}  // namespace StaticCondensationUtility
  
}  // namespace Kratos.
