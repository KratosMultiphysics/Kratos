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
            ElementType& rTheElement,
            MatrixType& rLeftHandSideMatrix,
            const std::vector<int> & rDofList)
        {
            KRATOS_TRY;
            const SizeType num_dofs_condensed = rDofList.size();
            const SizeType num_dofs_remaining = GetNumDofsElement(rTheElement) - num_dofs_condensed;
            const double numerical_limit = std::numeric_limits<double>::epsilon();

            std::vector<MatrixType> sub_matrices = CalculateSchurComplements(rTheElement, rLeftHandSideMatrix,rDofList);
            //1.) inverse K22
            MatrixType K_temp = ZeroMatrix(sub_matrices[3].size1());
            double detK22 = 0.00;
            MathUtils<double>::InvertMatrix(sub_matrices[3],K_temp,detK22);
            KRATOS_ERROR_IF(std::abs(detK22) < numerical_limit) << "Element " << rTheElement.Id()
                                                                << " is singular !" << std::endl;

            //2.) K_cond = K11 - K12*inv(K22)*K21
            K_temp = prod(K_temp,sub_matrices[2]);
            K_temp = prod(sub_matrices[1],K_temp);
            K_temp = sub_matrices[0]-K_temp;

            //3.) Fill rLeftHandSide to maintain same matrix size
            const std::vector<int> remaining_dof_list = CreateRemainingDofList(rTheElement, rDofList);
            rLeftHandSideMatrix.clear();

            SizeType dofA = 0;
            SizeType dofB = 0;
            for (SizeType i=0; i<num_dofs_remaining; ++i)
            {
                dofA = remaining_dof_list[i];
                for (SizeType j=0; j<num_dofs_remaining; ++j)
                {
                    dofB = remaining_dof_list[j];
                    rLeftHandSideMatrix(dofA,dofB) = K_temp(i,j);
                }
            }
            KRATOS_CATCH("")
        }

        std::vector<MatrixType> CalculateSchurComplements(
            ElementType& rTheElement,
            const MatrixType& rLeftHandSideMatrix,
            const std::vector<int> & rDofList)
        {
            KRATOS_TRY;
            // K11(0) K12(1)
            // K21(2) K22(3)		K22->dofs to be cond.
            // rDofList -> List of dofs to be condensed
            const std::vector<int> remaining_dof_list = CreateRemainingDofList(rTheElement, rDofList);
            const SizeType num_dofs_condensed = rDofList.size();
            const SizeType num_dofs_remaining = GetNumDofsElement(rTheElement)-num_dofs_condensed;

            KRATOS_ERROR_IF(num_dofs_remaining != remaining_dof_list.size()) << "unequal remaining dof size" << std::endl;

            std::vector<MatrixType> sub_matrices(4);
            sub_matrices[0] = ZeroMatrix(num_dofs_remaining, num_dofs_remaining);
            sub_matrices[1] = ZeroMatrix(num_dofs_remaining, num_dofs_condensed);
            sub_matrices[2] = ZeroMatrix(num_dofs_condensed, num_dofs_remaining);
            sub_matrices[3] = ZeroMatrix(num_dofs_condensed, num_dofs_condensed);

            FillSchurComplements(sub_matrices[0], rLeftHandSideMatrix, remaining_dof_list,
                                    remaining_dof_list, num_dofs_remaining, num_dofs_remaining);

            FillSchurComplements(sub_matrices[1], rLeftHandSideMatrix, remaining_dof_list,
                                    rDofList, num_dofs_remaining, num_dofs_condensed);

            FillSchurComplements(sub_matrices[2], rLeftHandSideMatrix, rDofList,
                                    remaining_dof_list, num_dofs_condensed, num_dofs_remaining);

            FillSchurComplements(sub_matrices[3], rLeftHandSideMatrix, rDofList,
                                    rDofList, num_dofs_condensed, num_dofs_condensed);

            return sub_matrices;
            KRATOS_CATCH("")
        }

        std::vector<int> CreateRemainingDofList(
            ElementType& rTheElement,
            const std::vector<int> & rDofList)
        {
            KRATOS_TRY;
            const SizeType num_dofs_condensed = rDofList.size();

            //fill remaining dofs
            std::vector<int> remaining_dofs_vec(0);
            for (SizeType i=0; i<GetNumDofsElement(rTheElement); ++i)
            {
                int current_dof = i;
                bool check = false;
                for (SizeType j = 0; j<num_dofs_condensed; ++j)
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

            for (SizeType i=0; i<rSizeA; ++i)
            {
                current_dof_a = rVecA[i];
                for (SizeType j=0; j<rSizeB; ++j)
                {
                    current_dof_b = rVecB[j];
                    Submatrix(i,j) = rLeftHandSideMatrix(current_dof_a, current_dof_b);
                }
            }
            KRATOS_CATCH("")
        }

        void ConvertingCondensation(
            ElementType& rTheElement,
            Vector& rLocalizedDofVector,
            Vector& rValues,
            const std::vector<int>& rDofList,
            const MatrixType& rLeftHandSideMatrix)
        {
            KRATOS_TRY;
            const double numerical_limit = std::numeric_limits<double>::epsilon();
            const std::vector<int> remaining_dof_list = CreateRemainingDofList(rTheElement, rDofList);
            const SizeType num_dofs_condensed = rDofList.size();

            const SizeType num_dofs_element = GetNumDofsElement(rTheElement);

            const SizeType num_dofs_remaining = num_dofs_element - num_dofs_condensed;
            std::vector<MatrixType> sub_matrices = CalculateSchurComplements(rTheElement, rLeftHandSideMatrix,rDofList);

            //1.) create u1
            Vector remaining_dofs_disp = ZeroVector(num_dofs_remaining);
            // Vector all_dofs_disp = ZeroVector(num_dofs_element);
            // rTheElement.GetValuesVector(all_dofs_disp);
            // rTheElement.LocalizeVector(all_dofs_disp); // localize global displacement -> element lvl
            // Note: "rLocalizedDofVector" is what was "all_dofs_disp" previously
            for (SizeType i=0; i<num_dofs_remaining; ++i) remaining_dofs_disp[i] = rLocalizedDofVector[remaining_dof_list[i]];

            //2.) inverse K22
            MatrixType K22_inv = ZeroMatrix(sub_matrices[3].size1());
            double detK22 = 0.00;
            MathUtils<double>::InvertMatrix(sub_matrices[3], K22_inv, detK22);

            KRATOS_ERROR_IF(std::abs(detK22) < numerical_limit) << "Element " << rTheElement.Id() << " is singular !" << std::endl;

            //3.) u2=inv(K22)*(F2-K21*u1),F2=0->u2=-inv(K22)*K21*u1
            Vector CondensedDofsDisp = ZeroVector(num_dofs_condensed);
            CondensedDofsDisp = prod(sub_matrices[2],remaining_dofs_disp);
            CondensedDofsDisp = -prod(K22_inv,CondensedDofsDisp);

            //4.) Fill rValues to maintain same matrix size
            rValues = ZeroVector(num_dofs_element);
            for (int i=0; i<static_cast<int>(num_dofs_element); ++i)
            {
                bool check = false;
                //check if dof i is condensed
                for (SizeType j=0;j<num_dofs_condensed;++j)
                {
                    if (i == rDofList[j])
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
                    for (SizeType j=0; j<num_dofs_remaining; ++j)
                    {
                        if (i == remaining_dof_list[j])
                        {
                            rValues[i] = remaining_dofs_disp[j];
                            break;
                        }
                    }
                }
            }
            // rTheElement.GlobalizeVector(rValues); // globalize local displacements -> global lvl
            KRATOS_CATCH("")
        }

        SizeType GetNumDofsElement(ElementType& rTheElement)
        {
            Vector all_dof_values = Vector(0);
            rTheElement.GetValuesVector(all_dof_values);
            return all_dof_values.size();
        }

	}  // namespace StaticCondensationUtility

}  // namespace Kratos.
