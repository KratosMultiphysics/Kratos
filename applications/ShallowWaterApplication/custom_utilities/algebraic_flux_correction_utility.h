//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


#ifndef KRATOS_ALGEBRAIC_FLUX_CORRECTION_UTILITY_H_INCLUDED
#define KRATOS_ALGEBRAIC_FLUX_CORRECTION_UTILITY_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
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

/// Utility to estimate the time step in terms of the courant number.
/** The velocity can be the sum of the convective velocity and the wave speed
*/
class KRATOS_API(SHALLOW_WATER_APPLICATION) AlgebraicFluxCorrectionUtility
{
public:
    ///@name Type Definitions
    ///@{

    using size_t = std::size_t;

    typedef DenseVector<Vector> ElementalValuesVectorType;

    typedef std::vector<Dof<double>::Pointer> DofsVectorType;

    typedef std::vector<DofsVectorType> ElementalDofsVectorType;

    /// Pointer definition of AlgebraicFluxCorrectionUtility
    KRATOS_CLASS_POINTER_DEFINITION(AlgebraicFluxCorrectionUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    AlgebraicFluxCorrectionUtility(ModelPart& rThisModelPart, Parameters ThisParameters);

    /// Destructor.
    ~AlgebraicFluxCorrectionUtility(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    int Check();

    void InitializeCorrection();

    void GetHighOrderValues();

    void GetLowOrderValues();

    void ApplyCorrection();

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
    std::string Info() const
    {
        return "AlgebraicFluxCorrectionUtility";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
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

    ModelPart& mrModelPart;
    int mRebuildLevel;
    std::vector<const Variable<double>*> mLimitingVariables;
    ElementalValuesVectorType mHighOrderValues;
    ElementalValuesVectorType mLowOrderValues;
    ElementalValuesVectorType mPreviousValues;
    ElementalValuesVectorType mAlgebraicFluxCorrections;
    ElementalValuesVectorType mElementalMassMatrices;
    ElementalDofsVectorType mElementalDofs;
    Vector mElementalLimiters;
    std::vector<Vector> mLimitingNodalHighOrderValues;
    std::vector<Vector> mLimitingNodalLowOrderValues;
    ElementalValuesVectorType mDofsLowOrderValues;

    const double mEpsilon = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    const Parameters GetDefaultParameters() const;

    void GenerateVariablesList(const std::vector<std::string>& rVariablesNames);

    void ResizeNodalAndElementalVectors();

    void InitializeNonhistoricalVariables();

    void AssembleElementalMassMatrices();

    void GetElementalDofList();

    void ComputeElementalAlgebraicFluxCorrections();

    void ComputeLimiters();

    void ComputeLimiterSingleVariable(
        const Variable<double>& rVariable,
        const Vector& rHighOrderValues,
        const Vector& rLowOrderValues);

    void AssembleLimitedCorrections();

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    AlgebraicFluxCorrectionUtility& operator=(AlgebraicFluxCorrectionUtility const& rOther);

    /// Copy constructor.
    AlgebraicFluxCorrectionUtility(AlgebraicFluxCorrectionUtility const& rOther);


    ///@}

}; // Class AlgebraicFluxCorrectionUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                AlgebraicFluxCorrectionUtility& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const AlgebraicFluxCorrectionUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ALGEBRAIC_FLUX_CORRECTION_UTILITY_H_INCLUDED  defined
