//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
//

#pragma once

// Project includes
#include "processes/process.h"


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

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

/// Utility to calculate the nodal gradient separately for the positive and negative sides of the zero level-set function (interface)

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) CalculateLevelsetConsistentNodalGradientProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LumpedInterfacePositiveNegativePressureGradient
    KRATOS_CLASS_POINTER_DEFINITION(CalculateLevelsetConsistentNodalGradientProcess);

    /// Auxiliary container to be used as TLS
    typedef std::tuple<BoundedMatrix<double,3,2>, array_1d<double,3>, array_1d<double,3>, array_1d<double,3>, array_1d<double,3>, array_1d<double,3>> TLSContainerType2D;
    typedef std::tuple<BoundedMatrix<double,4,3>, array_1d<double,4>, array_1d<double,4>, array_1d<double,4>, array_1d<double,3>, array_1d<double,4>> TLSContainerType3D;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor with separate parameters
     *
     * @param rModelPart Complete model part (including boundaries) for the process to operate on
     */
    CalculateLevelsetConsistentNodalGradientProcess(
        ModelPart& rModelPart);

    /// Constructor with Kratos parameters.
    CalculateLevelsetConsistentNodalGradientProcess(
        ModelPart& rModelPart,
        Parameters Parameters);

    /// Constructor with Kratos model
    CalculateLevelsetConsistentNodalGradientProcess(
        Model& rModel,
        Parameters Parameters);

    /// Destructor.
    ~CalculateLevelsetConsistentNodalGradientProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Execution of the process
     */
    void Execute() override;

    /**
     * @brief This method provides the default parameters
     */
    const Parameters GetDefaultParameters() const override;

    // ///@}
    // ///@name Inquiry
    // ///@{

    // ///@}
    // ///@name Input and output
    // ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "CalculateLevelsetConsistentNodalGradientProcess";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "CalculateLevelsetConsistentNodalGradientProcess";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}

    ///@}
    ///@name Friends
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    // Reference to the model part
    ModelPart& mrModelPart;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

        bool IsSplit(const Vector& rDistances);

        TLSContainerType2D SetTLSContainer2D();

        TLSContainerType3D SetTLSContainer3D();

        std::function<void(Element& rElement, TLSContainerType2D& rTLSContainer)> GetScalarNodalGradientElementFunction2D();

        std::function<void(Element& rElement, TLSContainerType3D& rTLSContainer)> GetScalarNodalGradientElementFunction3D();

        template<class TTLSContainer>
        void CalculateScalarNodalGradientElementContribution(
            Element& rElement,
            TTLSContainer& rTLSContainer);

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

}; // Class CalculateLevelsetConsistentNodalGradientProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.
