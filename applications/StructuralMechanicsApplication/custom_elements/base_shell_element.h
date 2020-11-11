//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Based on the work of Massimo Petracca and Peter Wilson
//

#if !defined(KRATOS_BASE_SHELL_ELEMENT_H_INCLUDED)
#define KRATOS_BASE_SHELL_ELEMENT_H_INCLUDED


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

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

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

    using SizeType = std::size_t;

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
    ~BaseShellElement() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
    * ELEMENTS inherited from this class have to implement next
    * Create and Clone methods: MANDATORY
    */

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

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
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

    IntegrationMethod mIntegrationMethod = GeometryData::GI_GAUSS_2;

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

    void BaseInitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo);

    void BaseFinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo);

    void BaseInitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo);

    void BaseFinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo);

    virtual void SetupOrientationAngles();

    void CheckVariables() const;
    void CheckDofs() const;
    void CheckProperties(const ProcessInfo& rCurrentProcessInfo) const;
    void CheckSpecificProperties() const;

    /**
    * computes the local axis of the element (for visualization)
    * @param rVariable: the variable to select the output
    * @param rOutput: the computed local axis
    * @param rpCoordinateTransformation: the coordinate-transformation to be used for computing the local axis
    */
    template <typename T>
    void ComputeLocalAxis(const Variable<array_1d<double, 3> >& rVariable,
                          std::vector<array_1d<double, 3> >& rOutput,
                          const T& rpCoordinateTransformation) const
    {
        const SizeType num_gps = GetNumberOfGPs();
        if (rOutput.size() != num_gps) {
            rOutput.resize(num_gps);
        }

        for (IndexType i=1; i<num_gps; ++i) {
            noalias(rOutput[i]) = ZeroVector(3);
        }

        const auto localCoordinateSystem(rpCoordinateTransformation->CreateLocalCoordinateSystem());
        if (rVariable == LOCAL_AXIS_1) {
            noalias(rOutput[0]) = localCoordinateSystem.Vx();
        } else if (rVariable == LOCAL_AXIS_2) {
            noalias(rOutput[0]) = localCoordinateSystem.Vy();
        } else if (rVariable == LOCAL_AXIS_3) {
            noalias(rOutput[0]) = localCoordinateSystem.Vz();
        } else {
            KRATOS_ERROR << "Wrong variable: " << rVariable.Name() << "!" << std::endl;
        }
    }

    /**
    * computes the local material axis of the element (for visualization)
    * @param rVariable: the variable to select the output
    * @param rOutput: the computed local material axis
    * @param rpCoordinateTransformation: the coordinate-transformation to be used for computing the local material axis
    */
    template <typename T>
    void ComputeLocalMaterialAxis(const Variable<array_1d<double, 3> >& rVariable,
                                  std::vector<array_1d<double, 3> >& rOutput,
                                  const T& rpCoordinateTransformation) const
    {
        const double mat_angle = Has(MATERIAL_ORIENTATION_ANGLE) ? GetValue(MATERIAL_ORIENTATION_ANGLE) : 0.0;

        const SizeType num_gps = GetNumberOfGPs();
        if (rOutput.size() != num_gps) {
            rOutput.resize(num_gps);
        }

        for (IndexType i=1; i<num_gps; ++i) {
            noalias(rOutput[i]) = ZeroVector(3);
        }

        const auto localCoordinateSystem(rpCoordinateTransformation->CreateLocalCoordinateSystem());

        const auto eZ = localCoordinateSystem.Vz();

        if (rVariable == LOCAL_MATERIAL_AXIS_1) {
            const auto q = QuaternionType::FromAxisAngle(eZ(0), eZ(1), eZ(2), mat_angle);
            q.RotateVector3(localCoordinateSystem.Vx(), rOutput[0]);
        } else if (rVariable == LOCAL_MATERIAL_AXIS_2) {
            const auto q = QuaternionType::FromAxisAngle(eZ(0), eZ(1), eZ(2), mat_angle);
            q.RotateVector3(localCoordinateSystem.Vy(), rOutput[0]);
        } else if (rVariable == LOCAL_MATERIAL_AXIS_3) {
            noalias(rOutput[0]) = eZ;
        } else {
            KRATOS_ERROR << "Wrong variable: " << rVariable.Name() << "!" << std::endl;
        }
    }


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

#endif // KRATOS_BASE_SHELL_ELEMENT_H_INCLUDED  defined
