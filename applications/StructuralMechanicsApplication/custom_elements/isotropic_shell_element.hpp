// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Massimo Petracca
//

#if !defined(KRATOS_ISOTROPIC_SHELL_ELEMENT_H_INCLUDED )
#define  KRATOS_ISOTROPIC_SHELL_ELEMENT_H_INCLUDED



// System includes


// External includes

// Project includes
#include "includes/element.h"

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

/// Short class definition.
/** Detail class definition.
*/
class IsotropicShellElement
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of IsotropicShellElement
    KRATOS_CLASS_POINTER_DEFINITION( IsotropicShellElement );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IsotropicShellElement(IndexType NewId, GeometryType::Pointer pGeometry);
    IsotropicShellElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    ~IsotropicShellElement() override;


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param pGeom The pointer to the geometry of the element
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo) override;

    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo) override;

    void GetValuesVector(Vector& values, int Step) override;
    void GetFirstDerivativesVector(Vector& values, int Step = 0) override;
    void GetSecondDerivativesVector(Vector& values, int Step = 0) override;

    void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<double >& rVariable, std::vector<double>& Output, const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo) override;

    void Initialize() override;

    void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override;

    int  Check(const ProcessInfo& rCurrentProcessInfo) override;

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
    array_1d< BoundedMatrix<double,3,3> , 3 > mTs;
    BoundedMatrix<double,3,3> mTE0;

    array_1d< array_1d<double,3>, 3> rot_oldit;

    double mOrientationAngle;


    ///@}
    ///@name Private Operators
    ///@{
    void CalculateLocalGlobalTransformation(
        double& x12,
        double& x23,
        double& x31,
        double& y12,
        double& y23,
        double& y31,
        array_1d<double,3>& v1,
        array_1d<double,3>& v2,
        array_1d<double,3>& v3,
        double& area
    );

    void CalculateMembraneB(
        BoundedMatrix<double,9,3>& B,
        const double&  beta0,
        const double& loc1,
        const double& loc2,
        const double& loc3,
        const double& x12,
        const double& x23,
        const double& x31,
        const double& y12,
        const double& y23,
        const double& y31
    );


    void CalculateBendingB(
        BoundedMatrix<double,9,3>& Bb,
        const double& loc2,
        const double& loc3,
        const double& x12,
        const double& x23,
        const double& x31,
        const double& y12,
        const double& y23,
        const double& y31
    );

    void CalculateMembraneContribution(
        const BoundedMatrix<double,9,3>& Bm,
        const BoundedMatrix<double,3,3>& Em,
        BoundedMatrix<double,9,9>& Km
    );


    void AssembleMembraneContribution(
        const BoundedMatrix<double,9,9>& Km,
        const double& coeff,
        BoundedMatrix<double,18,18>& Kloc_system
    );

    void CalculateBendingContribution(
        const BoundedMatrix<double,9,3>& Bb,
        const BoundedMatrix<double,3,3>& Eb,
        BoundedMatrix<double,9,9>& Kb
    );

    void AssembleBendingContribution(
        const BoundedMatrix<double,9,9>& Kb,
        const double& coeff,
        BoundedMatrix<double,18,18>& Kloc_system
    );

    void CalculateGaussPointContribution(
        BoundedMatrix<double,18,18>& Kloc_system ,
        const BoundedMatrix<double,3,3>& Em,
        const BoundedMatrix<double,3,3>& Eb,
        const double& weight,
        const double& h, /*thickness*/
        const double& loc1, /*local coords*/
        const double& loc2,
        const double& loc3,
        const double& x12,
        const double& x23,
        const double& x31,
        const double& y12,
        const double& y23,
        const double& y31
    );

    double CalculateBeta(
        const BoundedMatrix<double,3,3>& Em
    );

    void CalculateMembraneElasticityTensor(
        BoundedMatrix<double,3,3>& Em,
        const double& h
    );

    void CalculateBendingElasticityTensor(
        BoundedMatrix<double,3,3>& Eb,
        const double& h );

    void CalculateAllMatrices(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
    );

    Vector& CalculateVolumeForce(
        Vector& rVolumeForce
    );

    void AddBodyForce(
        const double& h,
        const double& Area,
        const Vector& VolumeForce,
        VectorType& rRightHandSideVector
    );

    void RotateToGlobal(
        const array_1d<double,3>& v1,
        const array_1d<double,3>& v2,
        const array_1d<double,3>& v3,
        const BoundedMatrix<double,18,18>& Kloc_system,
        Matrix& rLeftHandSideMatrix
    );

    void RotateToGlobal(
        const array_1d<double,3>& v1,
        const array_1d<double,3>& v2,
        const array_1d<double,3>& v3,
        const BoundedMatrix<double,18,18>& Kloc_system,
        Matrix& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector
    );



    void NicePrint(const Matrix& A);

    void  AddVoigtTensorComponents(
        const double local_component,
        array_1d<double,6>& v,
        const array_1d<double,3>& a,
        const array_1d<double,3>& b
    );

    void CalculateAndAddKg(
        MatrixType& LHS,
        BoundedMatrix<double,18,18>& rWorkMatrix,
        const double& x12,
        const double& x23,
        const double& x31,
        const double& y12,
        const double& y23,
        const double& y31,
        const array_1d<double,3>& v1,
        const array_1d<double,3>& v2,
        const array_1d<double,3>& v3,
        const double& A
    );

    void CalculateKg_GaussPointContribution(
        BoundedMatrix<double,18,18>& Kloc_system ,
        const BoundedMatrix<double,3,3>& Em,
        const double& weight,
        const double& h, /*thickness*/
        const double& loc1, /*local coords*/
        const double& loc2,
        const double& loc3,
        const double& x12,
        const double& x23,
        const double& x31,
        const double& y12,
        const double& y23,
        const double& y31,
        const array_1d<double,9>& membrane_disp
    );

    void CalculateLocalShapeDerivatives(
        double alpha,
        BoundedMatrix<double,2,9>& DNu_loc ,
        BoundedMatrix<double,2,9>& DNv_loc ,
        BoundedMatrix<double,2,9>& DNw_loc ,
        const double& a, /*local coords*/ //loc1
        const double& b, //loc2
        const double& c, //loc3
        const double& x12,
        const double& x23,
        const double& x31,
        const double& y12,
        const double& y23,
        const double& y31
    );

    void CalculateProjectionOperator(
        BoundedMatrix<double,18,18>& rProjOperator,
        const double& x12,
        const double& x23,
        const double& x31,
        const double& y12,
        const double& y23,
        const double& y31
    );

    void ApplyProjection(
        BoundedMatrix<double,18,18>& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        BoundedMatrix<double,18,18>& rWorkMatrix,
        array_1d<double,18>& rWorkArray,
        const BoundedMatrix<double,18,18>& rProjOperator
    );

    void UpdateNodalReferenceSystem(
        const double& x12,
        const double& x23,
        const double& x31,
        const double& y12,
        const double& y23,
        const double& y31
    );

    void SaveOriginalReference(
        const array_1d<double,3>& v1,
        const array_1d<double,3>& v2,
        const array_1d<double,3>& v3
    );

    void CalculatePureDisplacement(
        Vector& values,
        const array_1d<double,3>& v1,
        const array_1d<double,3>& v2,
        const array_1d<double,3>& v3
    );

    void CalculatePureMembraneDisplacement(
        array_1d<double,9>& values,
        const array_1d<double,3>& v1,
        const array_1d<double,3>& v2,
        const array_1d<double,3>& v3
    );

    void CalculatePureBendingDisplacement(
        array_1d<double,9>& values,
        const array_1d<double,3>& v1,
        const array_1d<double,3>& v2,
        const array_1d<double,3>& v3
    );

    void InvertMatrix(
        const BoundedMatrix<double,3,3>& InputMatrix,
        BoundedMatrix<double,3,3>& InvertedMatrix,
        double& InputMatrixDet
    );

    void SetupOrientationAngles();


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization
    IsotropicShellElement() {}

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer,  Element )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer,  Element )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //IsotropicShellElement& operator=(const IsotropicShellElement& rOther);

    /// Copy constructor.
    //IsotropicShellElement(const IsotropicShellElement& rOther);


    ///@}

}; // Class IsotropicShellElement

/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				IsotropicShellElement& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				const IsotropicShellElement& rThis)
{
rThis.PrintInfo(rOStream);
rOStream << std::endl;
rThis.PrintData(rOStream);

return rOStream;
}*/
///@}

}  // namespace Kratos.
#endif // KRATOS_ISOTROPIC_SHELL_ELEMENT_H_INCLUDED  defined


