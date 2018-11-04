//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//                   Lorenzo Gracia
//


#if !defined(KRATOS_SMALL_DISPLACEMENT_INTERFACE_ELEMENT_H_INCLUDED )
#define  KRATOS_SMALL_DISPLACEMENT_INTERFACE_ELEMENT_H_INCLUDED

// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"

// Application includes
#include "custom_utilities/poro_element_utilities.hpp"
#include "custom_utilities/interface_element_utilities.hpp"

#include "dam_application_variables.h"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class SmallDisplacementInterfaceElement : public Element
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( SmallDisplacementInterfaceElement );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    SmallDisplacementInterfaceElement() : Element() {}

    // Constructor 1
    SmallDisplacementInterfaceElement(IndexType NewId, GeometryType::Pointer pGeometry) : Element( NewId, pGeometry ) {}

    // Constructor 2
    SmallDisplacementInterfaceElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : Element( NewId, pGeometry, pProperties )
    {
        // Lobatto integration method with the integration points located at the "mid plane nodes" of the interface
        mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
    }

    // Destructor
    virtual ~SmallDisplacementInterfaceElement() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    void Initialize() override;

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo ) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo ) override;

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetValuesVector(Vector& rValues, int Step = 0) override;

    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) override;

    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void SetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable, std::vector<ConstitutiveLaw::Pointer>& rValues,const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct ElementVariables
    {
        //Properties variable
        double Density;

        //Nodal variables;
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        array_1d<double,TNumNodes*TDim> VolumeAcceleration;

        //General elemental variables
        BoundedMatrix<double,TDim, TDim> RotationMatrix;
        array_1d<double,TDim> VoigtVector;

        //Variables computed at each GP

        //Constitutive Law parameters
        Vector StrainVector;
        Vector StressVector;
        Matrix ConstitutiveMatrix;
        Vector Np;
        Matrix GradNpT;
        Matrix F;
        double detF;
        //Auxiliary Variables
        BoundedMatrix<double,TDim, TNumNodes*TDim> Nu;
        array_1d<double,TDim> BodyAcceleration;
        double IntegrationCoefficient;
        double JointWidth;
        BoundedMatrix<double,TNumNodes*TDim,TNumNodes*TDim> UMatrix;
        BoundedMatrix<double,TDim,TDim> DimMatrix;
        BoundedMatrix<double,TNumNodes*TDim,TDim> UDimMatrix;
        array_1d<double,TNumNodes*TDim> UVector;
    };

    // Member Variables
    GeometryData::IntegrationMethod mThisIntegrationMethod;
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;
    std::vector<double> mInitialGap;
    std::vector<bool> mIsOpen;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateInitialGap(const GeometryType& Geom);

    void ExtrapolateGPJointWidth (const std::vector<double>& JointWidthContainer);

    void CalculateStiffnessMatrix( MatrixType& rStiffnessMatrix, const ProcessInfo& CurrentProcessInfo );

    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo );

    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo );

    void InitializeElementVariables(ElementVariables& rVariables,ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                    const GeometryType& Geom, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo);

    void CalculateRotationMatrix(BoundedMatrix<double,TDim,TDim>& rRotationMatrix, const GeometryType& Geom);

    void CalculateJointWidth(double& rJointWidth,const double& NormalRelDisp,const double& MinimumJointWidth,const unsigned int& GPoint);

    void CheckAndCalculateJointWidth(double& rJointWidth,ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                    double& rNormalRelDisp,const double& MinimumJointWidth,const unsigned int& GPoint);

    void CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& weight, const double& detJ);


    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);


    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, ElementVariables& rVariables);


    void CalculateOutputDoubles( std::vector<double>& rOutput, const std::vector<double>& GPValues );

    template< class TValueType >
    void CalculateOutputValues( std::vector<TValueType>& rOutput, const std::vector<TValueType>& GPValues );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    }


}; // Class SmallDisplacementInterfaceElement

} // namespace Kratos

#endif // KRATOS_SMALL_DISPLACEMENT_INTERFACE_ELEMENT_H_INCLUDED  defined
