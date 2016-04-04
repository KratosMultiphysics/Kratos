//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_U_PW_SMALL_STRAIN_ELEMENT_H_INCLUDED )
#define  KRATOS_U_PW_SMALL_STRAIN_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/U_Pw_element.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPwSmallStrainElement : public UPwElement<TDim,TNumNodes>
{

public:
    
    KRATOS_CLASS_POINTER_DEFINITION( UPwSmallStrainElement );

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPwElement<TDim,TNumNodes>::mThisIntegrationMethod;
    using UPwElement<TDim,TNumNodes>::mConstitutiveLawVector;
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    UPwSmallStrainElement(IndexType NewId = 0) : UPwElement<TDim,TNumNodes>( NewId ) {}
    
    /// Constructor using an array of nodes
    UPwSmallStrainElement(IndexType NewId, const NodesArrayType& ThisNodes) : UPwElement<TDim,TNumNodes>(NewId, ThisNodes) {}
    
    /// Constructor using Geometry
    UPwSmallStrainElement(IndexType NewId, GeometryType::Pointer pGeometry) : UPwElement<TDim,TNumNodes>(NewId, pGeometry) {}

    /// Constructor using Properties
    UPwSmallStrainElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : UPwElement<TDim,TNumNodes>( NewId, pGeometry, pProperties ) {}

    /// Copy Constructor
    UPwSmallStrainElement(UPwSmallStrainElement const& rOther) : UPwElement<TDim,TNumNodes>(rOther) {}

    /// Destructor
    virtual ~UPwSmallStrainElement() {}

    /// Assignment operator.
    UPwSmallStrainElement & operator=(UPwSmallStrainElement const& rOther)
    {
        UPwElement<TDim,TNumNodes>::operator=(rOther);
        
        return *this;
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;
    
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const;
    
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;
    
    int Check(const ProcessInfo& rCurrentProcessInfo);
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo);
    
    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rOutput, const ProcessInfo& rCurrentProcessInfo);
    
    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo);
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
        
    struct ElementVariables
    {
        ///Properties variables
        double DynamicViscosityInverse;
        double FluidDensity;
        double Density;
        double BiotCoefficient;
        double BiotModulusInverse;
        boost::numeric::ublas::bounded_matrix<double,TDim, TDim> PermeabilityMatrix;
        
        ///ProcessInfo variables
        double NewmarkCoefficientU;
        double NewmarkCoefficientP;
    
        ///Nodal variables
        array_1d<double,TNumNodes> PressureVector;
        array_1d<double,TNumNodes> DtPressureVector;
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        array_1d<double,TNumNodes*TDim> VelocityVector;
        array_1d<double,TNumNodes*TDim> VolumeAcceleration;

        ///General elemental variables
        Vector VoigtVector;
        
        ///Variables computed at each GP
        Matrix B;
        ///Constitutive Law parameters
        Vector StrainVector;
        Vector StressVector;
        Matrix ConstitutiveMatrix;
        Vector Np;
        Matrix GradNpT;
        Matrix F;
        double detF;
        ///Auxiliary Variables
        boost::numeric::ublas::bounded_matrix<double,TDim, TNumNodes*TDim> Nu;
        array_1d<double,TDim> BodyAcceleration;
        double IntegrationCoefficient;
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*TDim,TNumNodes*TDim> UMatrix;
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*TDim,TNumNodes> UPMatrix;
        boost::numeric::ublas::bounded_matrix<double,TNumNodes,TNumNodes*TDim> PUMatrix;
        boost::numeric::ublas::bounded_matrix<double,TNumNodes,TNumNodes> PMatrix;
        Matrix UVoigtMatrix;
        boost::numeric::ublas::bounded_matrix<double,TNumNodes,TDim> PDimMatrix;
        array_1d<double,TNumNodes*TDim> UVector;
        array_1d<double,TNumNodes> PVector;
    };
    
    /// Member Variables
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    
    void CaculateStiffnessMatrix( MatrixType& rStiffnessMatrix, const ProcessInfo& CurrentProcessInfo );
    
    
    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo );

    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo );

    void InitializeElementVariables(ElementVariables& rVariables,ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                    const GeometryType& Geom, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo);
    
    void CalculateBMatrix(Matrix& rB, const Matrix& GradNpT);

    
    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);
    
    void CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);
    
    void CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);
    
    void CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);
    
    void CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    
    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables);
    
    void CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector, ElementVariables& rVariables);
    
    void CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, ElementVariables& rVariables);
    
    void CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector, ElementVariables& rVariables);
    
    void CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables);
    
    void CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables);
    
    void CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    /// Serialization
    
    friend class Serializer;
    
    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    }


}; // Class UPwSmallStrainElement

} // namespace Kratos

#endif // KRATOS_U_PW_SMALL_STRAIN_ELEMENT_H_INCLUDED  defined
