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
#include "includes/process_info.h"
#include "includes/constitutive_law.h"

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
 * @class ComputeCauchyStressHistoryUtility
 * @ingroup StructuralMechanicsApplication
 * @brief Node Search
 * @details This class provides several methods to perform paralelized loops in c++ to compute stresses and store them in a matrix
 * within a given radius.
 * @author Manuel Messmer
 */

class ComputeCauchyStressHistoryUtility
{
    public:
    ///@name Type Definitions
    ///@{

    /// Node type definition
    using NodeType = Node;

    /// Geometry definitions
    using GeometryType = Geometry<NodeType>;

    /// Pointer definition of ComputeCauchyStressHistoryUtility
    KRATOS_CLASS_POINTER_DEFINITION(ComputeCauchyStressHistoryUtility);


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ComputeCauchyStressHistoryUtility() {}

    /// Destructor.
    ~ComputeCauchyStressHistoryUtility(){}

    /*
    This method computes the stress history of a CL from a strain one. The CL needs a geometry (element provided)
    in softening type CLs.
    */
    Matrix ComputeStressHistory(
        ConstitutiveLaw::Pointer pConstitutiveLaw,
        const Matrix& rStrainHistory,
        const GeometryType& rGeometry,
        const Properties& rProperties//,
        //const ConstitutiveLaw::StressMeasure& rStressMeasure = ConstitutiveLaw::StressMeasure::StressMeasure_Cauchy
        )
    {
        const auto strain_size = pConstitutiveLaw->GetStrainSize();
        const auto dim = pConstitutiveLaw->WorkingSpaceDimension();
        const auto n_nodes = rGeometry.size();
        const ConstitutiveLaw::StressMeasure &r_stress_measure = ConstitutiveLaw::StressMeasure::StressMeasure_Cauchy;

        // Aux shape functions and derivatives
        Vector N(n_nodes);
        N.clear();
        Matrix DN_DX(n_nodes, dim);
        DN_DX.clear();

        // Unused F in small strains...
        Matrix F(dim, dim);
        F.clear();
        for (int i = 0; i<dim; ++i)
            F(i, i) = 1.0;
        double detF = 1.0;

        // Constitutive tangent matrix of one step
        Matrix C(strain_size, strain_size);
        C.clear();
        auto process_info = ProcessInfo();

        // The stress vector of one step
        Vector stress_vector(strain_size);
        stress_vector.clear();

        // The strain vector of one step
        Vector strain_vector(strain_size);
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

        // Check the mat props
        pConstitutiveLaw->Check(rProperties, rGeometry, process_info);

        // We initialize the material
        pConstitutiveLaw->ResetMaterial(rProperties, rGeometry, N);
        pConstitutiveLaw->InitializeMaterial(rProperties, rGeometry, N);

        const unsigned int n_steps = rStrainHistory.size1();

        const bool requires_initialize_resp = pConstitutiveLaw->RequiresInitializeMaterialResponse();
        const bool requires_finalize_resp   = pConstitutiveLaw->RequiresFinalizeMaterialResponse();

        for (unsigned int step = 0; step < n_steps; ++step) {
            noalias(strain_vector) = row(rStrainHistory, step);

            if (requires_initialize_resp)
                pConstitutiveLaw->InitializeMaterialResponse(cl_parameters, r_stress_measure);

            pConstitutiveLaw->CalculateMaterialResponse(cl_parameters, r_stress_measure);
            noalias(row(stress_history, step)) = cl_parameters.GetStressVector();

            if (requires_finalize_resp)
                pConstitutiveLaw->FinalizeMaterialResponse(cl_parameters, r_stress_measure);
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

    }; // Class ComputeCauchyStressHistoryUtility

///@}

///@} addtogroup block

}  // namespace Kratos.