//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_U_PW_SMALL_STRAIN_INTERFACE_ELEMENT_H_INCLUDED )
#define  KRATOS_U_PW_SMALL_STRAIN_INTERFACE_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/U_Pw_element.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/interface_element_utilities.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPwSmallStrainInterfaceElement : public UPwElement<TDim,TNumNodes>
{

public:
    
    KRATOS_CLASS_POINTER_DEFINITION( UPwSmallStrainInterfaceElement );

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
    UPwSmallStrainInterfaceElement(IndexType NewId = 0) : UPwElement<TDim,TNumNodes>( NewId ) {}
    
    /// Constructor using an array of nodes
    UPwSmallStrainInterfaceElement(IndexType NewId, const NodesArrayType& ThisNodes) : UPwElement<TDim,TNumNodes>(NewId, ThisNodes) {}
    
    /// Constructor using Geometry
    UPwSmallStrainInterfaceElement(IndexType NewId, GeometryType::Pointer pGeometry) : UPwElement<TDim,TNumNodes>(NewId, pGeometry) {}

    /// Constructor using Properties
    UPwSmallStrainInterfaceElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : UPwElement<TDim,TNumNodes>( NewId, pGeometry, pProperties ) 
    {
        /// Lobatto integration method with the integration points located at the "mid plane nodes" of the interface
        mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
    }

    /// Copy Constructor
    UPwSmallStrainInterfaceElement(UPwSmallStrainInterfaceElement const& rOther) : UPwElement<TDim,TNumNodes>(rOther), mInitialGap(rOther.mInitialGap) {}

    /// Destructor
    virtual ~UPwSmallStrainInterfaceElement() {}

    /// Assignment operator.
    UPwSmallStrainInterfaceElement & operator=(UPwSmallStrainInterfaceElement const& rOther)
    {
        UPwElement<TDim,TNumNodes>::operator=(rOther);
        
        mInitialGap.resize(rOther.mInitialGap.size(),false);
        
        for(unsigned int i = 0; i<mInitialGap.size(); i++)
            mInitialGap[i] = rOther.mInitialGap[i];
        
        return *this;
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;
    
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const;
    
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;
    
    int Check(const ProcessInfo& rCurrentProcessInfo);
    
    void Initialize();
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
    
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);
    
    void GetValueOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rValues, const ProcessInfo& rCurrentProcessInfo);
    
    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);
    
    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rOutput, const ProcessInfo& rCurrentProcessInfo);
    
    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo);
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct SFGradAuxVariables
    {
        array_1d<double,TDim> GlobalCoordinatesGradients;
        array_1d<double,TDim> LocalCoordinatesGradients;
        
        boost::numeric::ublas::bounded_matrix<double,TNumNodes,TDim-1> ShapeFunctionsNaturalGradientsMatrix;
        boost::numeric::ublas::bounded_matrix<double,TDim-1,TDim-1> LocalCoordinatesGradientsMatrix;
        boost::numeric::ublas::bounded_matrix<double,TDim-1,TDim-1> LocalCoordinatesGradientsInvMatrix;
        boost::numeric::ublas::bounded_matrix<double,TNumNodes,TDim-1> ShapeFunctionsGradientsMatrix;
    };
        
    struct InterfaceElementVariables
    {        
        ///Properties variables
        double DynamicViscosityInverse;
        double FluidDensity;
        double Density;
        double BiotCoefficient;
        double BiotModulusInverse;
        
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
        boost::numeric::ublas::bounded_matrix<double,TDim, TDim> RotationMatrix;
        array_1d<double,TDim> VoigtVector;
        
        ///Variables computed at each GP
        
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
        boost::numeric::ublas::bounded_matrix<double,TDim, TDim> LocalPermeabilityMatrix;
        array_1d<double,TDim> BodyAcceleration;
        double IntegrationCoefficient;
        double JointWidth;
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*TDim,TNumNodes*TDim> UMatrix;
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*TDim,TNumNodes> UPMatrix;
        boost::numeric::ublas::bounded_matrix<double,TNumNodes,TNumNodes*TDim> PUMatrix;
        boost::numeric::ublas::bounded_matrix<double,TNumNodes,TNumNodes> PMatrix;
        boost::numeric::ublas::bounded_matrix<double,TDim,TDim> DimMatrix;
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*TDim,TDim> UDimMatrix;
        boost::numeric::ublas::bounded_matrix<double,TNumNodes,TDim> PDimMatrix;
        array_1d<double,TNumNodes*TDim> UVector;
        array_1d<double,TNumNodes> PVector;
    };
    
    /// Member Variables
    
    Vector mInitialGap;
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateInitialGap(const GeometryType& Geom);
    
    
    void CaculateStiffnessMatrix( MatrixType& rStiffnessMatrix, const ProcessInfo& CurrentProcessInfo );
    
    
    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo );

    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo );

    void InitializeElementVariables(InterfaceElementVariables& rVariables,ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                    const GeometryType& Geom, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo);
                                    
    void CalculateRotationMatrix(boost::numeric::ublas::bounded_matrix<double,TDim,TDim>& rRotationMatrix, const GeometryType& Geom);
    
    void CalculateJointWidth(double& rJointWidth,const double& NormalRelDisp,const double& MinimumJointWidth,const unsigned int& GPoint);
    
    void CheckAndCalculateJointWidth(double& rJointWidth,ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                    double& rNormalRelDisp,const double& MinimumJointWidth,const unsigned int& GPoint);

    template< class TMatrixType >
    void CalculateShapeFunctionsGradients(TMatrixType& rGradNpT, SFGradAuxVariables& rAuxVariables,const Matrix& Jacobian, 
                                            const boost::numeric::ublas::bounded_matrix<double,TDim,TDim>& RotationMatrix,
                                            const Matrix& DN_De,const Matrix& Ncontainer, const double& JointWidth,const unsigned int& GPoint);
                                            
                                    
    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, InterfaceElementVariables& rVariables);
    
    void CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, InterfaceElementVariables& rVariables);
    
    void CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix, InterfaceElementVariables& rVariables);
    
    void CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, InterfaceElementVariables& rVariables);
    
    void CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix, InterfaceElementVariables& rVariables);

    
    void CalculateAndAddRHS(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables);
    
    void CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables);
    
    void CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables);
    
    void CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables);
    
    void CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables);
    
    void CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables);
    
    void CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables);


    void InterpolateOutputDoubles( std::vector<double>& rOutput, const std::vector<double>& GPValues );
    
    template< class TValueType >
    void InterpolateOutputValues( std::vector<TValueType>& rOutput, const std::vector<TValueType>& GPValues );
    
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


}; // Class UPwSmallStrainInterfaceElement

} // namespace Kratos

#endif // KRATOS_U_PW_SMALL_STRAIN_INTERFACE_ELEMENT_H_INCLUDED  defined 
