//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h"
#include "utilities/coordinate_transformation_utilities.h"

namespace Kratos
{

///@name Kratos Globals
///@{


///@}
///@name Type Definitions
///@{

///@}


///@name  Enum's
///@{


///@}
///@name  Functions
///@{



///@}
///@name Kratos Classes
///@{

/// Scheme for the solution of problems involving a slip condition.
/** This Scheme is a reimplementation of ResidualBasedIncrementalUpdateStaticScheme that can be used to
  * apply slip conditions along the edges of the model. The problem is solved using rotated coordinate axes on
  * the nodes where this condition will be applied, with the first coordinate of the rotated system aligned with the
  * normal to the boundary on each of these nodes.
  */
template<class TSparseSpace, class TDenseSpace>
class IncrementalUpdateRotationScheme : public ResidualBasedIncrementalUpdateStaticSchemeSlip<TSparseSpace,TDenseSpace>
{

public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( IncrementalUpdateRotationScheme);

    using BaseSchemeType = Scheme<TSparseSpace,TDenseSpace>;

    using BaseType = ResidualBasedIncrementalUpdateStaticSchemeSlip<TSparseSpace,TDenseSpace>;

    using TDataType = typename BaseType::TDataType;

    using DofsArrayType = typename BaseType::DofsArrayType;

    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    using TSystemVectorType = typename BaseType::TSystemVectorType;

    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;

    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    using RotationToolType = CoordinateTransformationUtils<LocalSystemMatrixType,LocalSystemVectorType,double>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit IncrementalUpdateRotationScheme() : BaseType()
    {
    }

    /**
     * @brief Constructor. The pseudo static scheme (parameters)
     * @param ThisParameters Configuration parameters
     */
    explicit IncrementalUpdateRotationScheme(Parameters ThisParameters)
        : BaseType()
    {
        // Check if the position of the rotation DOF is given
        KRATOS_WARNING_IF("IncrementalUpdateRotationScheme", !ThisParameters.Has("rotation_dof_position"))
            << "'rotation_dof_position' is not provided. Assuming it to be the first DOF in the node." << std::endl;

        // Validate default parameters
        Parameters default_parameters = Parameters(R"({
            "name" : "IncrementalUpdateRotationScheme",
            "block_size" : 0,
            "domain_size" : 0,
            "rotation_dof_position" : 0,
            "rotation_flag_name" : "SLIP"
        })");
        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        // Get rotation data
        const int block_size = ThisParameters["block_size"].GetInt();
        const int domain_size = ThisParameters["domain_size"].GetInt();
        const int rotation_dof_position = ThisParameters["rotation_dof_position"].GetInt();
        const auto& r_flag = KratosComponents<Flags>::Get(ThisParameters["rotation_flag_name"].GetString());

        // Check the provided parameters
        KRATOS_ERROR_IF(block_size == 0) << "'block_size' needs to be provided." << std::endl;
        KRATOS_ERROR_IF(block_size == 0) << "'domain_size' needs to be provided." << std::endl;

        // Set the pointer to the rotation tool
        BaseType::mpRotationTool = Kratos::make_shared<RotationToolType>(domain_size, block_size, rotation_dof_position, r_flag);
    }

    /// Destructor.
    ~IncrementalUpdateRotationScheme() override
    {
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    typename BaseSchemeType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<IncrementalUpdateRotationScheme<TSparseSpace, TDenseSpace>>(ThisParameters);
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "incremental_update_rotation_scheme";
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "IncrementalUpdateRotationScheme";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

}; // class

///@}
///@name Type Definitions
///@{

///@}

///@} // group

}  // namespace Kratos

