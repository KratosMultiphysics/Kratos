//
//   Project Name:        Kratos GeoMechanics
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 25 Jan 2017 $
//   Revision:            $Revision: 1.6 $
//
//


#if !defined(KRATOS_PARTIALLY_SATURATED_SOILS_KINEMATIC_LINEAR_H_INCLUDED )
#define  KRATOS_PARTIALLY_SATURATED_SOILS_KINEMATIC_LINEAR_H_INCLUDED



// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "custom_phase_laws/univariate_phase_law.h"


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
/**
Partially-saturated soils 3-phase element
REF:
+   Nagel's thesis
+   TR1
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPaPwSmallStrainElement : public Element
{

public:
    ///@name Type Definitions
    ///@{

    typedef Element BaseType;

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef ConstitutiveLaw ConstitutiveLawType;

    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    typedef Element::GeometryType GeometryType;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPaPwSmallStrainElement );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    UPaPwSmallStrainElement( IndexType NewId, GeometryType::Pointer pGeometry );
    UPaPwSmallStrainElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Destructor.
    virtual ~UPaPwSmallStrainElement();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /****************************CORE COMPUTATION********************************************************/

    Element::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const override;

    Element::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo ) const final;

    void GetDofList( DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo ) const final;

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) final;

    void CalculateDampingMatrix(MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo) final;

    void FinalizeSolutionStep( const ProcessInfo& CurrentProcessInfo ) final;

    int Check( const ProcessInfo& rCurrentProcessInfo ) const final;

    /**************************POST PROCESSING**********************************************************/

    void GetValuesVector( Vector& values, int Step ) const final;

    void CalculateOnIntegrationPoints( const Variable<double>& rVariable,
            std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo ) final;

    void CalculateOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable,
            std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo ) final;

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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "UPaPwSmallStrainElement #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "UPaPwSmallStrainElement #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;


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

    virtual void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                       const ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag,
                       bool CalculateResidualVectorFlag );

    //************************************************************************************
    //************************************************************************************
    //SATURATION AND ITS DERIVATIVES
    //************************************************************************************
    //************************************************************************************

    virtual double GetSaturation( const double& capillaryPressure ) const;

    virtual double GetDerivativeDSaturationDpc( const double& capillaryPressure ) const;

    virtual double GetSecondDerivativeD2SaturationDpc2( const double& capillaryPressure ) const;

    //************************************************************************************
    //************************************************************************************
    //PRIMARY VARIABLES AND THEIR DERIVATIVES
    //************************************************************************************
    //************************************************************************************

    virtual double GetDerivativeDWaterPressureDt( const Vector& N_PRESS ) const;

    virtual double GetDerivativeDAirPressureDt( const Vector& N_PRESS ) const;

    virtual double GetWaterPressure( const Vector& N ) const;

    virtual double GetAirPressure( const Vector& N ) const;

    virtual Vector GetGradientWater( const Matrix& DN_DX_PRESS ) const;

    virtual Vector GetGradientAir( const Matrix& DN_DX_PRESS ) const;

    //************************************************************************************
    //************************************************************************************
    //POROSITY AND ITS DERIVATIVES
    //************************************************************************************
    //************************************************************************************

    double GetPorosity( const Matrix& DN_DX_DISP ) const;

    double GetDerivativeDPorosityDDivU( const Matrix& DN_DX_DISP ) const;

    //************************************************************************************
    //************************************************************************************
    //DISPLACEMENT AND ITS DERIVATIVES
    //************************************************************************************
    //************************************************************************************

    double GetDivU( const Matrix& DN_DX_DISP ) const;

    double GetDerivativeDDivUDt( const Matrix& DN_DX_DISP ) const;

    //************************************************************************************
    //************************************************************************************
    //STRESSES, STRAINS AND CONSTITUTIVE MODEL
    //************************************************************************************
    //************************************************************************************

    void CalculateBoperator( Matrix& B_Operator, const Matrix& DN_DX ) const;

    void CalculateStrain( const Matrix& B, const Matrix& Displacements, Vector& StrainVector ) const;

    //************************************************************************************
    //************************************************************************************
    //FLUID FLOW AND ITS DERIVATIVES
    //************************************************************************************
    //************************************************************************************

    virtual void GetWaterValuesOnIntegrationPoints( const std::size_t& ipoint, double& permeability_water, double& density_water ) const;

    Vector GetFlowWater( const Vector& grad_water,
            const double& saturation, const double& permeability_water, const double& density_water ) const;

    Vector GetDerivativeDWaterFlowDpw( const Vector& flow_water,
            const double& saturation, const double& DS_Dpc ) const;

    Vector GetDerivativeDWaterFlowDpa( const Vector& flow_water,
            const double& saturation, const double& DS_Dpc ) const;

    double GetDerivativeDWaterFlowDGradpw( const double& saturation, const double& permeability_water, const double& density_water ) const;

    void GetAirValuesOnIntegrationPoints( const std::size_t& ipoint, double& permeability_air ) const;

    Vector GetFlowAir( const Vector& grad_air, const double& saturation, const double& permeability_air, const double& density_air ) const;

    Vector GetDerivativeDAirFlowDpa( const Vector& flow_air, const double& saturation, const double& DS_Dpc, const double& permeability_air, const double& density_air, const double& bulk_air ) const;

    Vector GetDerivativeDAirFlowDpw( const Vector& flow_air, const double& saturation, const double& DS_Dpc ) const;

    double GetDerivativeDAirFlowDGradpa( const double& saturation, const double& permeability_air, const double& density_air ) const;

    //************************************************************************************
    //************************************************************************************
    //DENSITY
    //************************************************************************************
    //************************************************************************************

    double GetAveragedDensity( const double& saturation, const double& porosity,
        const double& density_soil, const double& density_water, const double& density_air ) const;

    double GetDensityAir( const double& airPressure ) const;

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{

    // A protected default constructor necessary for serialization
    UPaPwSmallStrainElement() {}

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    UnivariatePhaseLaw::Pointer mpSWCCLaw;
    UnivariatePhaseLaw::Pointer mpRelativePermeabilityWaterLaw;
    UnivariatePhaseLaw::Pointer mpRelativePermeabilityAirLaw;
    UnivariatePhaseLaw::Pointer mpGasLaw;

    GeometryType::Pointer mpSubGeometry;

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    IntegrationMethod mThisIntegrationMethod;

    Matrix mInitialDisp;
    std::vector<double> mReferencePressures;
    bool mIsStabilized;

    ///@}
    ///@name Private Operations
    ///@{

    void InitializeSubGeometry();

    void InitializeMaterial(const ProcessInfo& rCurrentProcessInfo);

    //************************************************************************************
    //************************************************************************************
    //STRESSES, STRAINS AND CONSTITUTIVE MODEL
    //************************************************************************************
    //************************************************************************************

    void CalculateEffectiveStress( Vector& StressVector, Matrix& tanC_W, Matrix& tanC_A,
            const double& waterPressure, const double& airPressure, const double& saturation, const double& DS_Dpc ) const;

    void CalculateStressAndTangentialStiffness( const int& PointNumber,
            const Vector& StrainVector, Vector& StressVector, Matrix& tanC_U,
            const double& waterPressure, const double& airPressure,
            const double& saturation, const double& DS_Dpc,
            Matrix& tanC_W, Matrix& tanC_A,
            const ProcessInfo& rCurrentProcessInfo );

    //************************************************************************************
    //************************************************************************************
    //CALCULATE FORCEVECTORS DISPLACEMENT
    //************************************************************************************
    //************************************************************************************

    void CalculateBodyForcesToRHSVectorU( Vector& R, const Vector& N_DISP, const double& density ) const;

    void CalculateInternalForcesToRHSU( Vector& R, const Matrix& B_Operator, const Vector& StressVector ) const;

    //************************************************************************************
    //************************************************************************************
    //CALCULATE STIFFNESS MATRICES DISPLACEMENT
    //************************************************************************************
    //************************************************************************************

    void CalculateStiffnessMatrixUU( Matrix& Help_K_UU, const Matrix& tanC, const Matrix& B_Operator,
                                    const Matrix& DN_DX_DISP, const Vector& N_DISP,
                                    const double& Drho_DdivU ) const;

    void CalculateStiffnessMatrixUW( Matrix& Help_K_UW, const Matrix& tanC_W, const Matrix& DN_DX_DISP,
                                    const Vector& N_DISP, const Vector& N_PRESS, const double& Drho_Dpw ) const;

    void CalculateStiffnessMatrixUA( Matrix& Help_K_UA, const Matrix& tanC_A, const Matrix& DN_DX_DISP,
                                    const Vector& N_DISP, const Vector& N_PRESS, const double& Drho_Dpa ) const;

    //************************************************************************************
    //************************************************************************************
    //CALCULATE FORCEVECTORS WATER
    //************************************************************************************
    //************************************************************************************

    void CalculateInternalForcesToRHSW( Vector& Help_R_U, const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, const Vector& N_PRESS,
                                        const double& saturation, const double& porosity,
                                        const double& DS_Dpc, const double& Dpc_Dt, const double& divU_Dt,
                                        const Vector& flow_water ) const;

    //************************************************************************************
    //************************************************************************************
    //CALCULATE STIFFNESS MATRICES WATER
    //************************************************************************************
    //************************************************************************************

    void CalculateStiffnessMatrixWU( Matrix& Help_K_WU, const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, const Vector& N_PRESS,
                                     const double& DS_Dpc, const double& Dpc_Dt, const double& Dn_DdivU ) const;

    void CalculateStiffnessMatrixWW( Matrix& Help_K_WW, const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, const Vector& N_PRESS,
                                     const double& porosity, const double& DS_Dpc, const double& Dpc_Dt, const double& D2S_Dpc2,
                                     const double& divU_Dt, const double& Dflow_waterDgradpw, const Vector& Dflow_waterDpw ) const;

    void CalculateStiffnessMatrixWA( Matrix& Help_K_WA, const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, const Vector& N_PRESS,
                                     const double& porosity, const double& DS_Dpc, const double& Dpc_Dt, const double& D2S_Dpc2,
                                     const double& divU_Dt, const Vector& Dflow_waterDpa ) const;

    //************************************************************************************
    //************************************************************************************
    //CALCULATE DAMPING MATRICES WATER
    //************************************************************************************
    //************************************************************************************

    void CalculateDampingMatrixWU( Matrix& Help_D_WU, const Matrix& DN_DX_DISP, const Vector& N_PRESS,
                                   const double& saturation ) const;

    void CalculateDampingMatrixWW( Matrix& Help_D_WW, const Matrix& DN_DX_DISP, const Vector& N_PRESS,
                                   const double& porosity, const double& DS_Dpc ) const;

    void CalculateDampingMatrixWA( Matrix& Help_D_WA, const Matrix& DN_DX_DISP, const Vector& N_PRESS,
                                   const double& porosity, const double& DS_Dpc ) const;

    //************************************************************************************
    //************************************************************************************
    //CALCULATE FORCEVECTORS WATER
    //************************************************************************************
    //************************************************************************************

    void CalculateInternalForcesToRHSA( Vector& Help_R_A, const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, const Vector& N_PRESS,
                                        const double& saturation, const double& porosity,  const double& density_air, const double& bulk_air,
                                        const double& DS_Dpc, const double& Dpc_Dt, const double& divU_Dt,
                                        const Vector& flow_air, const Vector& grad_air, const double& airPressure_Dt ) const;

    //************************************************************************************
    //************************************************************************************
    //CALCULATE STIFFNESS MATRICES AIR
    //************************************************************************************
    //************************************************************************************

    void CalculateStiffnessMatrixAU( Matrix& Help_K_AU, const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, const Vector& N_PRESS,
                                     const double& saturation, const double& DS_Dpc, const double& Dpc_Dt,
                                     const double& Dn_DdivU, const double& density_air, const double& bulk_air,
                                     const double& airPressure_Dt ) const;

    void CalculateStiffnessMatrixAW( Matrix& Help_K_AW, const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, const Vector& N_PRESS,
                                     const double& porosity, const double& DS_Dpc, const double& Dpc_Dt, const double& D2S_Dpc2,
                                     const double& divU_Dt, const double& density_air, const double& bulk_air,
                                     const double& airPressure_Dt, const Vector& Dflow_airDpw, const Vector& grad_air ) const;

    void CalculateStiffnessMatrixAA( Matrix& Help_K_AA, const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, const Vector& N_PRESS,
                                     const double& saturation,
                                     const double& porosity, const double& DS_Dpc, const double& Dpc_Dt, const double& D2S_Dpc2,
                                     const double& divU_Dt, const double& density_air, const double& bulk_air,
                                     const double& airPressure_Dt, const Vector& flow_air,
                                     const double& Dflow_airDgradpa, const Vector& Dflow_airDpa,
                                     const Vector& grad_air ) const;

    //************************************************************************************
    //************************************************************************************
    //CALCULATE DAMPING MATRICES AIR
    //************************************************************************************
    //************************************************************************************

    void CalculateDampingMatrixAU( Matrix& Help_D_AU, const Matrix& DN_DX_DISP, const Vector& N_PRESS,
                                   const double& saturation ) const;

    void CalculateDampingMatrixAW( Matrix& Help_D_AW, const Matrix& DN_DX_DISP, const Vector& N_PRESS,
                                   const double& porosity, const double& DS_Dpc ) const;

    void CalculateDampingMatrixAA( Matrix& Help_D_AA, const Matrix& DN_DX_DISP, const Vector& N_PRESS,
                                   const double& saturation, const double& porosity, const double& DS_Dpc,
                                   const double& density_air, const double& bulk_air ) const;

    //************************************************************************************
    //************************************************************************************
    //AUXILIARY FUNCTIONS
    //************************************************************************************
    //************************************************************************************

    void AssembleRHSFromSubVectors( VectorType& rRightHandSideVector,
            const unsigned int& dim,
            const unsigned int& num,
            const unsigned int& stride,
            const unsigned int& begin,
            const Vector& R_P ) const;

    void AssembleStiffnessFromSubMatrices( MatrixType& rLeftHandSideMatrix,
            const unsigned int& dim1,
            const unsigned int& num1,
            const unsigned int& stride1,
            const unsigned int& begin1,
            const unsigned int& dim2,
            const unsigned int& num2,
            const unsigned int& stride2,
            const unsigned int& begin2,
            const Matrix& K_PQ) const;

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer,  Element );
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //UPaPwSmallStrainElement& operator=(const UPaPwSmallStrainElement& rOther);

    /// Copy constructor.
    //UPaPwSmallStrainElement(const UPaPwSmallStrainElement& rOther);

    ///@}

}; // Class UPaPwSmallStrainElement

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
                   UPaPwSmallStrainElement& rThis);
*/

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const UPaPwSmallStrainElement& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
   return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_PARTIALLY_SATURATED_SOILS_KINEMATIC_LINEAR_H_INCLUDED defined


