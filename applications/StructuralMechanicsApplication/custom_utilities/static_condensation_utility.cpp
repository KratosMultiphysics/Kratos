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
#include "includes/element.h"
#include "static_condensation_utility.h"


namespace Kratos
{
	namespace StaticCondensationUtility
	{
		void CondenseLeftHandSide(
            const ElementType& rTheElement,
            MatrixType& rLeftHandSideMatrix,
			const std::vector<int> & rDofList)
        {
            KRATOS_TRY;
            const SizeType dofs_condensed = rDofList.size();
            const SizeType dofs_remaining = msElementSize-dofs_condensed;
            const double numerical_limit = std::numeric_limits<double>::epsilon();

            std::vector<MatrixType> SubMatrices = CalculateSchurComplements(rTheElement, rLeftHandSideMatrix,rDofList);
            //1.) inverse K22
            MatrixType K_temp = ZeroMatrix(SubMatrices[3].size1());
            double detK22 = 0.00;
            MathUtils<double>::InvertMatrix(SubMatrices[3],K_temp,detK22);
            KRATOS_ERROR_IF(std::abs(detK22) < numerical_limit) << "Element " 
                                                                << rTheElement.Id() << " is singular !" << std::endl;
            
            //2.) K_cond = K11 - K12*inv(K22)*K21 
            K_temp = prod(K_temp,SubMatrices[2]);
            K_temp = prod(SubMatrices[1],K_temp);
            K_temp = SubMatrices[0]-K_temp;

            //3.) Fill rLeftHandSide to maintain same matrix size
            const std::vector<int> RemainingDofList = CreateRemainingDofList(rTheElement, rDofList);
            rLeftHandSideMatrix.clear();

            SizeType dofA = 0;
            SizeType dofB = 0;
            for (SizeType i =0;i<dofs_remaining;++i)
            {
                dofA = RemainingDofList[i];
                for (SizeType j=0;j<dofs_remaining;++j)
                {
                    dofB = RemainingDofList[j];
                    rLeftHandSideMatrix(dofA,dofB) = K_temp(i,j);
                }
            }
            KRATOS_CATCH("")
        }

        std::vector<MatrixType> CalculateSchurComplements(
            const ElementType& rTheElement,
			const MatrixType& rLeftHandSideMatrix,
			const std::vector<int> & rDofList)
        {
            KRATOS_TRY;
            // K11(0) K12(1)	
            // K21(2) K22(3)		K22->dofs to be cond.
            // rDofList -> List of dofs to be condensed
            const std::vector<int> RemainingDofList = CreateRemainingDofList(rTheElement, rDofList);
            const SizeType dofs_condensed = rDofList.size();
            const SizeType dofs_remaining = msElementSize-dofs_condensed;

            KRATOS_ERROR_IF(dofs_remaining != RemainingDofList.size()) << "unequal remaining dof size" << std::endl;

            std::vector<MatrixType> SubMatrices(4);
            SubMatrices[0] = ZeroMatrix(dofs_remaining, dofs_remaining);
            SubMatrices[1] = ZeroMatrix(dofs_remaining, dofs_condensed);
            SubMatrices[2] = ZeroMatrix(dofs_condensed, dofs_remaining);
            SubMatrices[3] = ZeroMatrix(dofs_condensed, dofs_condensed);

            FillSchurComplements(SubMatrices[0],rLeftHandSideMatrix,RemainingDofList,RemainingDofList,dofs_remaining,dofs_remaining);
            FillSchurComplements(SubMatrices[1],rLeftHandSideMatrix,RemainingDofList,rDofList,dofs_remaining,dofs_condensed);
            FillSchurComplements(SubMatrices[2],rLeftHandSideMatrix,rDofList,RemainingDofList,dofs_condensed,dofs_remaining);
            FillSchurComplements(SubMatrices[3],rLeftHandSideMatrix,rDofList,rDofList,dofs_condensed,dofs_condensed);

            return SubMatrices;
            KRATOS_CATCH("") 
        }

		std::vector<int> CreateRemainingDofList(
            const ElementType& rTheElement,
			const std::vector<int> & rDofList)
        {
            KRATOS_TRY;
            const SizeType dofs_condensed = rDofList.size();

            //fill remaining dofs
            int current_dof = 0;
            bool check = false;
            std::vector<int> remaining_dofs_vec(0);
            for (SizeType i=0;i<msElementSize;++i)
            {
                current_dof = i;
                check = false;
                for (SizeType j = 0;j<dofs_condensed;++j)
                {
                    if (current_dof == rDofList[j]) check = true;
                }
                if (check) continue;
                else remaining_dofs_vec.push_back(current_dof);
            }
            return remaining_dofs_vec;

            KRATOS_CATCH("");
        }

		void FillSchurComplements(
			MatrixType& Submatrix,
			const MatrixType& rLeftHandSideMatrix,
			const std::vector<int>& rVecA,
			const std::vector<int>& rVecB,
			const SizeType& rSizeA,
			const SizeType& rSizeB) //maybe inline
        {
            KRATOS_TRY;
            SizeType current_dof_a = 0;
            SizeType current_dof_b = 0;

            for (SizeType i=0;i<rSizeA;++i)
            {
                current_dof_a = rVecA[i];
                for (SizeType j=0;j<rSizeB;++j)
                {
                    current_dof_b = rVecB[j];				
                    Submatrix(i,j) = rLeftHandSideMatrix(current_dof_a,current_dof_b);
                }
            }
            KRATOS_CATCH("")
        }

		void ConvertingCondensation(
            const ElementType& rTheElement,
			Vector& rValues,
			const std::vector<int>& rDofList,
			const MatrixType& rLeftHandSideMatrix)
        {
            KRATOS_TRY;
            const double numerical_limit = std::numeric_limits<double>::epsilon();
            const std::vector<int> RemainingDofList = CreateRemainingDofList(rTheElement, rDofList);
            const SizeType dofs_condensed = rDofList.size();
            const SizeType dofs_remaining = msElementSize-dofs_condensed;
            std::vector<MatrixType> SubMatrices = CalculateSchurComplements(rTheElement, rLeftHandSideMatrix,rDofList);

            //1.) create u1 
            Vector RemainingDofsDisp = ZeroVector(dofs_remaining);
            Vector AllDofsDisp = ZeroVector(msElementSize);
            rTheElement.GetValuesVector(AllDofsDisp);
            rTheElement.LocalizeVector(AllDofsDisp); // localize global displacement -> element lvl
            for (SizeType i=0;i<dofs_remaining;++i) RemainingDofsDisp[i] = AllDofsDisp[RemainingDofList[i]];

            //2.) inverse K22
            MatrixType K22_inv = ZeroMatrix(SubMatrices[3].size1());
            double detK22 = 0.00;
            MathUtils<double>::InvertMatrix(SubMatrices[3],K22_inv,detK22);

            KRATOS_ERROR_IF(std::abs(detK22) < numerical_limit) << "Element " << rTheElement.Id() << " is singular !" << std::endl;

            //3.) u2=inv(K22)*(F2-K21*u1),F2=0->u2=-inv(K22)*K21*u1
            Vector CondensedDofsDisp = ZeroVector(dofs_condensed);
            CondensedDofsDisp = prod(SubMatrices[2],RemainingDofsDisp);
            CondensedDofsDisp = -prod(K22_inv,CondensedDofsDisp);

            //4.) Fill rValues to maintain same matrix size
            rValues = ZeroVector(msElementSize);
            bool check;
            for (int i=0;i<msElementSize;++i)
            {
                check = false;
                //check if dof i is condensed
                for (SizeType j=0;j<dofs_condensed;++j)
                {
                    if (i==rDofList[j])
                    {
                        rValues[i] = CondensedDofsDisp[j];
                        check = true;
                        break;
                    } 
                }

                if (check) continue; // found respective dof -> search for next dof
                //check remaining dofs
                else
                {
                    for (SizeType j=0;j<dofs_remaining;++j)
                    {
                        if (i==RemainingDofList[j])
                        {
                            rValues[i] = RemainingDofsDisp[j];
                            break;
                        }
                    }
                }
            }
            rTheElement.GlobalizeVector(rValues); // globalize local displacements -> global lvl
            KRATOS_CATCH("")
        }

	}  // namespace StaticCondensationUtility
  
}  // namespace Kratos.
