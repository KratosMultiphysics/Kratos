// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//


#if !defined(KRATOS_TOTAL_LAGRANGIAN_H_INCLUDED )
#define  KRATOS_TOTAL_LAGRANGIAN_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "custom_elements/base_solid_element.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"


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

/// Total Lagrangian element for 2D and 3D geometries.

/**
 * Implements a total Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 2D and 3D
 */

class TotalLagrangian
    : public BaseSolidElement
{
public:
    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// Counted pointer of TotalLagrangian
    KRATOS_CLASS_POINTER_DEFINITION(TotalLagrangian);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TotalLagrangian(IndexType NewId, GeometryType::Pointer pGeometry);
    TotalLagrangian(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~TotalLagrangian();

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{
    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
    //TODO: ADD THE OTHER CREATE FUNCTION
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

    void Initialize();

    void ResetConstitutiveLaw();

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

    void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& Output, const ProcessInfo& rCurrentProcessInfo);

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& Output, const ProcessInfo& rCurrentProcessInfo);

    void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo& rCurrentProcessInfo);

    void SetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);


    void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);

    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo);

    //std::string Info() const;

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
    //      virtual String Info() const;

    /// Print information about this object.
    //      virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    //      virtual void PrintData(std::ostream& rOStream) const;
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
	TotalLagrangian() : BaseSolidElement()
    {
    }

    /**
     * Calculates the elemental contributions
     * \f$ K^e = w\,B^T\,D\,B \f$ and
     * \f$ r^e \f$
     */
    virtual void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo,
                              bool CalculateStiffnessMatrixFlag,
                              bool CalculateResidualVectorFlag);
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

    ///@}
    ///@name Private Operators
    ///@{

    void CalculateAndAddKm(
        MatrixType& K,
        Matrix& B,
        Matrix& D,
        double weight);

    /**
     * Calculation of the Geometric Stiffness Matrix. Kg = dB * S
     */
    void CalculateAndAddKg(
        MatrixType& K,
        Matrix& DN_DX,
        Vector& StressVector,
        double weight
    );

    void CalculateBodyForces(
        Vector& BodyForce,
        const ProcessInfo& CurrentProcessInfo
    );

    void InitializeVariables();

    virtual void InitializeMaterial();

    double CalculateIntegrationWeight
    (double GaussPointWeight,
     double DetJ0);

    void CalculateAndAdd_ExtForceContribution(
        const Vector& N,
        const ProcessInfo& CurrentProcessInfo,
        Vector& BodyForce,
        VectorType& mResidualVector,
        double weight
    );

    void CalculateStrain(const Matrix& C,
                         Vector& StrainVector);

    void CalculateB(Matrix& B,
                    Matrix& F,
                    Matrix& DN_DX,
                    unsigned int StrainSize);


    void Comprobate_State_Vector(Vector& Result);

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization



    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);



    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    /// Assignment operator.
    //TotalLagrangian& operator=(const TotalLagrangian& rOther);
    /// Copy constructor.
    //TotalLagrangian(const TotalLagrangian& rOther);
    ///@}

}; // Class TotalLagrangian

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_TOTAL_LAGRANGIAN_H_INCLUDED  defined 
