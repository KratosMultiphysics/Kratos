//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Guillermo Casas (gcasas@cimne.upc.edu)
//
//

#ifndef KRATOS_DERIVATIVE_RECOVERY_UTILITY_H
#define KRATOS_DERIVATIVE_RECOVERY_UTILITY_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"

// Application includes


namespace Kratos
{
///@addtogroup SwimmingDEMApplication
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

class KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecoveryUtility
{
public:
    ///@name Type Definitions
    ///@{
    typedef Variable<double> ScalarVariableType;
    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > ComponentVariableType;
    typedef Variable<array_1d<double, 3> > VectorVariableType;
    /// Pointer definition of DerivativeRecoveryUtility
    KRATOS_CLASS_POINTER_DEFINITION(DerivativeRecoveryUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Kratos parameters.
    DerivativeRecoveryUtility(
        ModelPart& rModelPart,
        Parameters rParameters);

    /// Constructor with Kratos model
    DerivativeRecoveryUtility(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    virtual ~DerivativeRecoveryUtility() = default;

    ///@}
    ///@name Operators
    ///@{

    virtual void Initialize(){};
    virtual void Recover();

    ///@}
    ///@name Operations
    ///@{

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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "DerivativeRecoveryUtility" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "DerivativeRecoveryUtility";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

bool mStoreFullGradient;
ModelPart& mrModelPart;

virtual void AddTimeDerivative(const ScalarVariableType& rScalarVariable, const ScalarVariableType& rTimeDerivativeVariable){}

virtual void AddTimeDerivative(const ComponentVariableType& rScalarComponent, const ScalarVariableType& rTimeDerivativeVariable){}

virtual void AddTimeDerivative(const VectorVariableType& rVectorVariable, const VectorVariableType& rTimeDerivativeVariable){}

virtual void CalculateGradient(const ScalarVariableType& rScalarVariable, const VectorVariableType& rGradientVariable){}

virtual void CalculateGradient(const ComponentVariableType& rScalarComponent, const VectorVariableType& rGradientVariable){}

virtual void CalculateGradient(const VectorVariableType& rVectorVariable,
                               const VectorVariableType& rComponent0GradientVariable,
                               const VectorVariableType& rComponent1GradientVariable,
                               const VectorVariableType& rComponent2GradientVariable){}

virtual void CalculateDivergence(const VectorVariableType rVectorVariable, const ScalarVariableType& rDivergenceVariable){}

virtual void CalculateLaplacian(const ScalarVariableType& rScalarVariable, const ScalarVariableType& rLaplacianVariable){}

virtual void CalculateLaplacian(const ComponentVariableType& rScalarComponent, const ScalarVariableType& rLaplacianVariable){}

virtual void CalculateLaplacian(const VectorVariableType& rVectorVariable, const VectorVariableType& rLaplacianVariable){}

virtual void CalculateMaterialDerivative(const ScalarVariableType& rScalarVariable, const ScalarVariableType& rMaterialDerivativeVariable){}

virtual void CalculateMaterialDerivative(const ComponentVariableType& rScalarComponent, const ScalarVariableType& rMaterialDerivativeVariable){}

virtual void CalculateMaterialDerivative(const VectorVariableType& rScalarComponent, const VectorVariableType& rMaterialDerivativeVariable){}

virtual void CalculateRotational(const VectorVariableType rVectorVariable, const VectorVariableType& rRotationalVariable){}

virtual void CheckDefaultsAndSettings(Parameters rParameters){};

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Protected Operators
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

    /// Default constructor.
    DerivativeRecoveryUtility() = delete;

    /// Assignment operator.
    DerivativeRecoveryUtility& operator=(DerivativeRecoveryUtility const& rOther) = delete;

    /// Copy constructor.
    DerivativeRecoveryUtility(DerivativeRecoveryUtility const& rOther) = delete;

    ///@}

}; // Class DerivativeRecoveryUtility

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_DERIVATIVE_RECOVERY_UTILITY_H
