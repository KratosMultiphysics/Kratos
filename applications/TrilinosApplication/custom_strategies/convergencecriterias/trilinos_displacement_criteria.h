//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#pragma once

//  System includes

//  External includes

//  Project includes
#include "solving_strategies/convergencecriterias/displacement_criteria.h"

namespace Kratos
{

///@addtogroup TrilinosApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class TrilinosDisplacementCriteria
 * @ingroup TrilinosApplication
 * @brief MPI version of the DisplacementCriteria.
 * @details This is a convergence criteria that considers the increment on the solution as criteria. The reactions from the RHS are not computed in the solution
 * @see DisplacementCriteria
 * @author Jordi Cotela
 */
template<class TSparseSpace,
         class TDenseSpace
         >
class TrilinosDisplacementCriteria
    : public DisplacementCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosDisplacementCriteria
    KRATOS_CLASS_POINTER_DEFINITION( TrilinosDisplacementCriteria );

    /// The definition of the base ConvergenceCriteria
    using BaseConvergenceCriteriaType = ConvergenceCriteria< TSparseSpace, TDenseSpace >;

    /// The definition of the base DisplacementCriteria
    using BaseType = DisplacementCriteria<TSparseSpace, TDenseSpace>;

    /// The definition of the class type
    using ClassType = TrilinosDisplacementCriteria<TSparseSpace, TDenseSpace>;

    /// The definition of the data type
    using TDataType = typename BaseType::TDataType;

    /// The definition of the DoF data type
    using DofType = typename Node::DofType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit TrilinosDisplacementCriteria()
        : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param ThisParameters The configuration parameters
     */
    explicit TrilinosDisplacementCriteria(Kratos::Parameters ThisParameters)
        : BaseType()
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * @brief Constructor 2 arguments
     * @param NewRatioTolerance The ratio tolerance for the convergence.
     * @param AlwaysConvergedNorm The absolute tolerance for the convergence.
     */
    explicit TrilinosDisplacementCriteria(
        TDataType NewRatioTolerance,
        TDataType AlwaysConvergedNorm)
        : BaseType(NewRatioTolerance, AlwaysConvergedNorm)
    {
    }

    /**
     * @brief Copy constructor
     * @param rOther The criteria to be copied
     */
    explicit TrilinosDisplacementCriteria( TrilinosDisplacementCriteria const& rOther )
      :BaseType(rOther)
    {
    }

    /**
     * @brief Destructor.
     */
    ~TrilinosDisplacementCriteria() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create method
     * @param ThisParameters The configuration parameters
     */
    typename BaseConvergenceCriteriaType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "trilinos_displacement_criteria";
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                            : "trilinos_displacement_criteria",
            "displacement_relative_tolerance" : 1.0e-4,
            "displacement_absolute_tolerance" : 1.0e-9
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "TrilinosDisplacementCriteria";
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
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);
    }

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
}; //  Class ClassName

///@}
///@name Type Definitions
///@{

///@}

}  //  namespace Kratos.
