// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//                   Based on the work of Massimo Petracca and Peter Wilson
//

#pragma once


// System includes

// External includes


// Project includes
#include "includes/element.h"
#include "utilities/quaternion.h"
#include "custom_utilities/shell_cross_section.hpp"
#include "structural_mechanics_application_variables.h"


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

enum class ShellKinematics
{
    LINEAR,
    NONLINEAR_COROTATIONAL
};

///@}
///@name Kratos Classes
///@{

template <class TCoordinateTransformation>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) BaseShellElement
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    typedef Element BaseType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of BaseShellElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(BaseShellElement);

    typedef std::vector< ShellCrossSection::Pointer > CrossSectionContainerType;

    typedef Quaternion<double> QuaternionType;

    using CoordinateTransformationPointerType = Kratos::unique_ptr<TCoordinateTransformation>;

    using SizeType = std::size_t;

    using Vector3Type = array_1d<double, 3>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Constructor using Geometry
    */
    BaseShellElement(IndexType NewId,
                     GeometryType::Pointer pGeometry);

    /**
    * Constructor using Properties
    */
    BaseShellElement(IndexType NewId,
                     GeometryType::Pointer pGeometry,
                     PropertiesType::Pointer pProperties);

    /**
    * Destructor
    */
    ~BaseShellElement() override = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
    * this determines the elemental equation ID vector for all elemental
    * DOFs
    * @param rResult: the elemental equation ID vector
    * @param rCurrentProcessInfo: the current process info instance
    */
    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const override;

    /**
    * determines the elemental list of DOFs
    * @param ElementalDofList: the list of DOFs
    * @param rCurrentProcessInfo: the current process info instance
    */
    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override;


    void GetValuesVector(Vector& rValues, int Step = 0) const override;

    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override;

    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;

    void ResetConstitutiveLaw() override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
                                      std::vector<array_1d<double, 3> >& rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
                                      std::vector<ConstitutiveLaw::Pointer>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    // Calculate functions
    void Calculate(const Variable<Matrix >& rVariable,
                   Matrix& Output,
                   const ProcessInfo& rCurrentProcessInfo) override;


    /**
    * This method provides the place to perform checks on the completeness of the input
    * and the compatibility with the problem options as well as the contitutive laws selected
    * It is designed to be called only once (or anyway, not often) typically at the beginning
    * of the calculations, so to verify that nothing is missing from the input
    * or that no common error is found.
    * @param rCurrentProcessInfo
    * this method is: MANDATORY
    */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    /**
    * returns the used integration method. In the general case this is the
    * default integration method of the used geometry. I an other integration
    * method is used the method has to be overwritten within the element
    * @return default integration method of the used Geometry
    */
    IntegrationMethod GetIntegrationMethod() const override
    {
        return mIntegrationMethod;
    }
    ///@}
    ///@name Access
    ///@{

    void SetCrossSectionsOnIntegrationPoints(std::vector< ShellCrossSection::Pointer >& crossSections);


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    const Parameters GetSpecifications() const override;

    /// Turn back information as a string.
    virtual std::string Info() const override;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override;

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

    IntegrationMethod mIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;

    CoordinateTransformationPointerType mpCoordinateTransformation = nullptr;

    CrossSectionContainerType mSections; /*!< Container for cross section associated to each integration point */

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{
    /**
    * Protected empty constructor
    */
    BaseShellElement() : Element()
    {
    }

    SizeType GetNumberOfDofs() const;

    SizeType GetNumberOfGPs() const;

    /**
     * @brief This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix The LHS
     * @param rRightHandSideVector The RHS
     * @param rCurrentProcessInfo The current process info instance
     * @param CalculateStiffnessMatrixFlag The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag The flag to set if compute the RHS
     */
    virtual void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    );

    void SetupOrientationAngles();

    void CheckDofs() const;
    void CheckProperties(const ProcessInfo& rCurrentProcessInfo) const;
    void CheckSpecificProperties() const;

    /**
    * computes the local axis of the element (for visualization)
    * @param rVariable: the variable to select the output
    * @param rOutput: the computed local axis
    * @param rpCoordinateTransformation: the coordinate-transformation to be used for computing the local axis
    */
    void ComputeLocalAxis(const Variable<array_1d<double, 3> >& rVariable,
                          std::vector<array_1d<double, 3> >& rOutput) const;

    /**
    * computes the local material axis of the element (for visualization)
    * @param rVariable: the variable to select the output
    * @param rOutput: the computed local material axis
    * @param rpCoordinateTransformation: the coordinate-transformation to be used for computing the local material axis
    */
    void ComputeLocalMaterialAxis(const Variable<array_1d<double, 3> >& rVariable,
                                  std::vector<array_1d<double, 3> >& rOutput) const;

    // check if this function is really necessary
    void DecimalCorrection(Vector& a);

    /**
    * Returns the behavior of this shell (thin/thick)
    * @return the shell behavior
    */
    virtual ShellCrossSection::SectionBehaviorType GetSectionBehavior() const;

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
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

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

}; // Class BaseShellElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
