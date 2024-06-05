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
// #include <string>
// #include <iostream>
// #include "includes/model_part.h"
#include "includes/process_info.h"

// Include kratos definitions

// Project includes

// Configures

// External includes

namespace Kratos {
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class ComputeCauchyStressHistory
 * @ingroup StructuralMechanicsApplication
 * @brief Node Search
 * @details This class provides several methods to perform paralelized loops in c++ to compute stresses and store them in a matrix
 * within a given radius.
 * @author Manuel Messmer
 */

class ComputeCauchyStressHistory
{
    public:
    ///@name Type Definitions
    ///@{

    /// Node type definition
    using NodeType = Node;

    /// Geometry definitions
    using GeometryType = Geometry<NodeType>;

    /// Pointer definition of ComputeCauchyStressHistory
    KRATOS_CLASS_POINTER_DEFINITION(ComputeCauchyStressHistory);


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ComputeCauchyStressHistory() {}

    /// Destructor.
    ~ComputeCauchyStressHistory(){}

    /*
    This method computes the stress history of a CL from a strain one. The CL needs a geometry (element provided)
    in softening type CLs.
    */
    Matrix ComputeStressHistory(
        ConstitutiveLaw::Pointer pConstitutiveLaw,
        const Matrix& rStrainHistory,
        const GeometryType& rGeometry,
        const Properties& rProperties
        )
    {
        // Aux shape functions and derivatives
        Vector N(3);
        N.clear();
        Matrix DN_DX(3, 2);
        DN_DX.clear();

        // Unused F in small strains...
        Matrix F(2, 2);
        F.clear();
        F(0, 0) = 1.0;
        F(1, 1) = 1.0;
        double detF = 1.0;

        // Constitutive tangent matrix of one step
        Matrix C(3, 3);
        C.clear();
        auto process_info = ProcessInfo();

        // The stress vector of one step
        Vector stress_vector(3);
        stress_vector.clear();

        // The strain vector of one step
        Vector strain_vector(3);
        strain_vector.clear();

        // Init the stress history matrix
        Matrix stress_history = rStrainHistory;
        stress_history.clear();


        auto cl_parameters = ConstitutiveLaw::Parameters();
        cl_parameters.SetConstitutiveMatrix(C);
        cl_parameters.SetStrainVector(strain_vector);
        cl_parameters.SetStressVector(stress_vector);
        cl_parameters.SetDeformationGradientF(F);
        cl_parameters.SetDeterminantF(detF);
        cl_parameters.SetShapeFunctionsValues(N);
        cl_parameters.SetShapeFunctionsDerivatives(DN_DX);
        cl_parameters.SetProcessInfo(process_info);
        cl_parameters.SetMaterialProperties(rProperties);
        cl_parameters.SetElementGeometry(rGeometry);

        // Set constitutive law flags:
        Flags& cl_options = cl_parameters.GetOptions();
        cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        // We initialize the material
        pConstitutiveLaw->InitializeMaterial(rProperties, rGeometry, N);

        const unsigned int n_steps = rStrainHistory.size1();

        for (unsigned int step = 0; step < n_steps; ++step) {
            noalias(strain_vector) = row(rStrainHistory, step);
            pConstitutiveLaw->CalculateMaterialResponseCauchy(cl_parameters);
            noalias(row(stress_history, step)) = cl_parameters.GetStressVector();
            pConstitutiveLaw->FinalizeMaterialResponseCauchy(cl_parameters);
        }

        return stress_history;
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

    }; // Class ComputeCauchyStressHistory

///@}

///@} addtogroup block

}  // namespace Kratos.