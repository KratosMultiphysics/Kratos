//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:        Lorenzo Gracia $
//   Date:                $Date:              March 2016 $
//   Revision:            $Revision:                 1.0 $
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
#include "custom_utilities/element_utilities.hpp"
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
    
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;
    
    void Initialize();
    
    int Check(const ProcessInfo& rCurrentProcessInfo);
    
    void GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
    
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);
    
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo );

    void CalculateRightHandSide(VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo );
    
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    void GetSecondDerivativesVector(Vector& rValues, int Step = 0);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void SetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);
    
    void GetValueOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rValues, const ProcessInfo& rCurrentProcessInfo);
    
    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rOutput, const ProcessInfo& rCurrentProcessInfo);

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
        boost::numeric::ublas::bounded_matrix<double,TDim, TDim> RotationMatrix;
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
        boost::numeric::ublas::bounded_matrix<double,TDim, TNumNodes*TDim> Nu;
        array_1d<double,TDim> BodyAcceleration;
        double IntegrationCoefficient;
        double JointWidth;
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*TDim,TNumNodes*TDim> UMatrix;
        boost::numeric::ublas::bounded_matrix<double,TDim,TDim> DimMatrix;
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*TDim,TDim> UDimMatrix; //TODO: YA NO ES NECESARIA
        array_1d<double,TNumNodes*TDim> UVector;
    };
    
    // Member Variables
    
    GeometryData::IntegrationMethod mThisIntegrationMethod;
    
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;
    
    Vector mInitialGap;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void CalculateInitialGap(const GeometryType& Geom);
    
    void CalculateStiffnessMatrix( MatrixType& rStiffnessMatrix, const ProcessInfo& CurrentProcessInfo );    
    
    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo );

    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo );

    void InitializeElementVariables(ElementVariables& rVariables,ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                    const GeometryType& Geom, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo);
    
    void CalculateRotationMatrix(boost::numeric::ublas::bounded_matrix<double,TDim,TDim>& rRotationMatrix, const GeometryType& Geom);
    
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
    
    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    }


}; // Class SmallDisplacementInterfaceElement

} // namespace Kratos

#endif // KRATOS_SMALL_DISPLACEMENT_INTERFACE_ELEMENT_H_INCLUDED  defined 
