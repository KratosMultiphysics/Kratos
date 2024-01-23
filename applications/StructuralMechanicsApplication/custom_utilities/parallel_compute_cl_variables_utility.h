

// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    A. Cornejo
//
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include "includes/model_part.h"

// Include kratos definitions

// Project includes

// Configures

// External includes

namespace Kratos
{
    ///@addtogroup StructuralMechanicsApplication
    ///@{

    ///@name Kratos Classes
    ///@{

    /**
     * @class ParallelComputeCLVariablesUtility
     * @ingroup StructuralMechanicsApplication
     * @brief Node Search
     * @details This class provides several methods to perform paralelized loops in c++ to compute stresses and store them in a matrix
     * within a given radius.
     * @author Manuel Messmer
     */

    class ParallelComputeCLVariablesUtility
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of ParallelComputeCLVariablesUtility
        KRATOS_CLASS_POINTER_DEFINITION(ParallelComputeCLVariablesUtility);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        ParallelComputeCLVariablesUtility() {}

        /// Destructor.
        ~ParallelComputeCLVariablesUtility() {}

        // Matrix ComputeStresses(ModelPart &rModelPart, const Matrix &rStrainList, std::vector<ConstitutiveLaw::Pointer> &rCLlist, const Flags &rOptions)
        Matrix ComputeStresses(ModelPart &rModelPart, const Matrix &rStrainList, std::vector<std::vector<ConstitutiveLaw::Pointer>> &rCLlist, const Flags &rOptions)
        {

            // Vector N(3);
            Vector N(4);

            N.clear();
            Matrix DN_DX(3, 2);
            DN_DX.clear();

            Matrix F(2, 2);
            F.clear();
            F(0, 0) = 1.0;
            F(1, 1) = 1.0;
            double detF = 1.0;

            Matrix C(3, 3);
            C.clear();
            auto &r_process_info = rModelPart.GetProcessInfo();

            Matrix stress_list = rStrainList;
            // Matrix stress_list(rModelPart.NumberOfElements(), rModelPart.GetElement(1).GetIntegrationPoints().size() * 3);

            stress_list.clear();
            // for (int i = 0; i < static_cast<int>(rCLlist.size()); ++i)Kratos::Element::IndexType
            const auto& r_elements_begin = rModelPart.ElementsBegin();

#pragma omp for
            for (int i = 0; i < static_cast<int>(rModelPart.NumberOfElements()); ++i)
            {
                // auto &r_element = rModelPart.GetElement(i+1);
                auto r_element = r_elements_begin+i;
                auto r_geometry = r_element->GetGeometry();
                auto &r_integration_points = r_geometry.IntegrationPoints();
                for (int j = 0; j < static_cast<int>(r_integration_points.size()); ++j)
                {

                // We initialize the element-wise info
                Vector stress_vector(3);
                stress_vector.clear();

                Vector strain_vector(3);
                strain_vector[0] = rStrainList(i, j*3);
                strain_vector[1] = rStrainList(i, j*3 + 1);
                strain_vector[2] = rStrainList(i, j*3 + 2);


                auto cl_parameters = ConstitutiveLaw::Parameters();
                cl_parameters.SetOptions(rOptions);
                cl_parameters.SetDeformationGradientF(F);
                cl_parameters.SetDeterminantF(detF);
                cl_parameters.SetStrainVector(strain_vector);
                cl_parameters.SetStressVector(stress_vector);
                cl_parameters.SetConstitutiveMatrix(C);
                cl_parameters.SetElementGeometry(r_geometry); // Use the geometry of the current element
                cl_parameters.SetShapeFunctionsValues(N);

                cl_parameters.SetShapeFunctionsDerivatives(DN_DX);
                cl_parameters.SetProcessInfo(r_process_info);
                cl_parameters.SetMaterialProperties(rModelPart.GetProperties(0));

                // strain_vector[0] = rStrainList(i, 0);
                // strain_vector[1] = rStrainList(i, 1);
                // strain_vector[2] = rStrainList(i, 2);
                // cl_parameters.SetElementGeometry(rModelPart.GetElement(i+1).GetGeometry()); // Here we assume that the ith CL corresponds to the ith element...
                // rCLlist[i]->CalculateMaterialResponseCauchy(cl_parameters); // Here I assume one IP
                rCLlist[i][j]->CalculateMaterialResponseCauchy(cl_parameters); // Now we use the CL of the current element and integration point
                // stress_list(i, 0) = stress_vector[0];
                // stress_list(i, 1) = stress_vector[1];
                // stress_list(i, 2) = stress_vector[2];
                stress_list(i, j*3) = stress_vector[0];
                stress_list(i, j*3 + 1) = stress_vector[1];
                stress_list(i, j*3 + 2) = stress_vector[2];
                // rCLlist[i]->FinalizeMaterialResponseCauchy(cl_parameters);
                rCLlist[i][j]->FinalizeMaterialResponseCauchy(cl_parameters);
                }
            }

            // block_for_each(rModelPart.Elements(), [&](Element& rElement) {
            //     // We initialize the element-wise info
            //     Vector stress_vector(3);
            //     stress_vector.clear();

            //     Vector strain_vector(3);

            //     auto cl_parameters = ConstitutiveLaw::Parameters();
            //     cl_parameters.SetOptions(rOptions);
            //     cl_parameters.SetDeformationGradientF(F);
            //     cl_parameters.SetDeterminantF(detF);
            //     cl_parameters.SetStrainVector(strain_vector);
            //     cl_parameters.SetStressVector(stress_vector);
            //     cl_parameters.SetConstitutiveMatrix(C);
            //     cl_parameters.SetShapeFunctionsValues(N);
            //     cl_parameters.SetShapeFunctionsDerivatives(DN_DX);
            //     cl_parameters.SetProcessInfo(r_process_info);
            //     cl_parameters.SetMaterialProperties(rModelPart.GetProperties(0));
            //     int elem_id = rElement.Id();

            //     std::vector<ConstitutiveLaw::Pointer> cl_list;
            //     rElement.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, cl_list, r_process_info);
            //     strain_vector[0] = rStrainList(elem_id, 0);
            //     strain_vector[1] = rStrainList(elem_id, 1);
            //     strain_vector[2] = rStrainList(elem_id, 2);
            //     cl_parameters.SetElementGeometry(rElement.GetGeometry());

            //     cl_list[0]->CalculateMaterialResponseCauchy(cl_parameters); // Here I assume one IP

            //     stress_list(elem_id, 0) = stress_vector[0];
            //     stress_list(elem_id, 1) = stress_vector[1];
            //     stress_list(elem_id, 2) = stress_vector[2];

            //     cl_list[0]->FinalizeMaterialResponseCauchy(cl_parameters);

            // });

            return stress_list;
        }

        ///@}
        ///@name Operations
        ///@{

        ///@}
        ///@name Input and output
        ///@{

        ///@}

    private:
        ///@name Member Variables
        ///@{

        ///@}
        ///@name Un accessible methods
        ///@{

        ///@}

    }; // Class ParallelComputeCLVariablesUtility

    ///@}

    ///@} addtogroup block

} // namespace Kratos.