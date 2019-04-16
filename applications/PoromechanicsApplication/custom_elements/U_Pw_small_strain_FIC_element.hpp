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
//


#if !defined(KRATOS_U_PW_SMALL_STRAIN_FIC_ELEMENT_H_INCLUDED )
#define  KRATOS_U_PW_SMALL_STRAIN_FIC_ELEMENT_H_INCLUDED

// System includes
#include <cmath>

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/U_Pw_element.hpp"
#include "custom_elements/U_Pw_small_strain_element.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPwSmallStrainFICElement : public UPwSmallStrainElement<TDim,TNumNodes>
{

public:
    
    KRATOS_CLASS_POINTER_DEFINITION( UPwSmallStrainFICElement );

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPwElement<TDim,TNumNodes>::mThisIntegrationMethod;
    using UPwElement<TDim,TNumNodes>::mConstitutiveLawVector;
    typedef typename UPwSmallStrainElement<TDim,TNumNodes>::ElementVariables ElementVariables;
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    UPwSmallStrainFICElement(IndexType NewId = 0) : UPwSmallStrainElement<TDim,TNumNodes>( NewId ) {}
    
    /// Constructor using an array of nodes
    UPwSmallStrainFICElement(IndexType NewId, const NodesArrayType& ThisNodes) : UPwSmallStrainElement<TDim,TNumNodes>(NewId, ThisNodes) {}
    
    /// Constructor using Geometry
    UPwSmallStrainFICElement(IndexType NewId, GeometryType::Pointer pGeometry) : UPwSmallStrainElement<TDim,TNumNodes>(NewId, pGeometry) {}

    /// Constructor using Properties
    UPwSmallStrainFICElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : UPwSmallStrainElement<TDim,TNumNodes>( NewId, pGeometry, pProperties ) {}

    /// Destructor
    ~UPwSmallStrainFICElement() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;
    
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;
        
    void Initialize() override;
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;
    
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
        
    struct FICElementVariables
    {
        /// Properties variables
        double ShearModulus;
        
        /// Nodal variables
        array_1d< array_1d<double,TDim*TNumNodes> , TNumNodes > NodalShapeFunctionsGradients;
        
        /// General elemental variables
        double ElementLength;
        Matrix VoigtMatrix;
        
        /// Variables computed at each GP
        BoundedMatrix<double,TDim,TDim*TNumNodes> StrainGradients;
        array_1d<Vector,TNumNodes> ShapeFunctionsSecondOrderGradients;
        array_1d< array_1d<double,TDim> , TDim > DtStressGradients;
        array_1d< std::vector< array_1d<double,TDim> > , TDim > ConstitutiveTensorGradients;
        ///Auxiliary variables
        array_1d<double,TDim> DimVector;
        Matrix DimVoigtMatrix;
        BoundedMatrix<double,TDim,TDim*TNumNodes> DimUMatrix;
    };
    
    /// Member Variables
    
    array_1d< std::vector< array_1d<double,TNumNodes> > , TDim > mNodalConstitutiveTensor;
    
    array_1d< array_1d<double,TNumNodes> , TDim > mNodalDtStress;
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    
    void SaveGPConstitutiveTensor(array_1d<Matrix,TDim>& rConstitutiveTensorContainer, const Matrix& ConstitutiveMatrix, const unsigned int& GPoint);
    
    void SaveGPDtStress(Matrix& rDtStressContainer, const Vector& StressVector, const unsigned int& GPoint);
    
    void ExtrapolateGPConstitutiveTensor(const array_1d<Matrix,TDim>& ConstitutiveTensorContainer);
    
    void ExtrapolateGPDtStress(const Matrix& DtStressContainer);
    
    
    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo ) override;

    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo ) override;

    void InitializeFICElementVariables(FICElementVariables& rFICVariables, const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer,
                                        const GeometryType& Geom,const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo);
    
    void ExtrapolateShapeFunctionsGradients(array_1d< array_1d<double,TDim*TNumNodes> , TNumNodes >& rNodalShapeFunctionsGradients,
                                                const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer);
    
    void CalculateElementLength(double& rElementLength, const GeometryType& Geom);
    
    void InitializeSecondOrderTerms(FICElementVariables& rFICVariables);
        
    void CalculateShapeFunctionsSecondOrderGradients(FICElementVariables& rFICVariables, ElementVariables& rVariables);
    
    
    void CalculateAndAddLHSStabilization(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables, FICElementVariables& rFICVariables);

    void CalculateAndAddStrainGradientMatrix(MatrixType& rLeftHandSideMatrix,ElementVariables& rVariables, FICElementVariables& rFICVariables);
    
    void CalculateAndAddDtStressGradientMatrix(MatrixType& rLeftHandSideMatrix,ElementVariables& rVariables, FICElementVariables& rFICVariables);
    
    void CalculateConstitutiveTensorGradients(FICElementVariables& rFICVariables, const ElementVariables& Variables);
    
    void CalculateAndAddPressureGradientMatrix(MatrixType& rLeftHandSideMatrix,ElementVariables& rVariables, FICElementVariables& rFICVariables);
    

    void CalculateAndAddRHSStabilization(VectorType& rRightHandSideVector, ElementVariables& rVariables, FICElementVariables& rFICVariables);

    void CalculateAndAddStrainGradientFlow(VectorType& rRightHandSideVector,ElementVariables& rVariables, FICElementVariables& rFICVariables);
    
    void CalculateAndAddDtStressGradientFlow(VectorType& rRightHandSideVector,ElementVariables& rVariables, FICElementVariables& rFICVariables);
    
    void CalculateDtStressGradients(FICElementVariables& rFICVariables, const ElementVariables& Variables);
    
    void CalculateAndAddPressureGradientFlow(VectorType& rRightHandSideVector,ElementVariables& rVariables, FICElementVariables& rFICVariables);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    /// Serialization
    
    friend class Serializer;
    
    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    }

    /// Assignment operator.
    UPwSmallStrainFICElement & operator=(UPwSmallStrainFICElement const& rOther);

    /// Copy constructor.
    UPwSmallStrainFICElement(UPwSmallStrainFICElement const& rOther);

}; // Class UPwSmallStrainFICElement

} // namespace Kratos

#endif // KRATOS_U_PW_SMALL_STRAIN_FIC_ELEMENT_H_INCLUDED  defined
