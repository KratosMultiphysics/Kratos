//
//   Project Name:        KratosDamApplication  $
//   Last modified by:    $Author: lgracia		$
//   Date:                $Date: December 2016  $
//   Revision:            $Revision: 1.0        $
//
//

#if !defined(KRATOS_WAVE_EQUATION_ELEMENT_H_INCLUDED )
#define  KRATOS_WAVE_EQUATION_ELEMENT_H_INCLUDED

// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"

#include "dam_application_variables.h" 

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class  WaveEquationElement : public Element
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( WaveEquationElement );
        
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    WaveEquationElement(IndexType NewId = 0) : Element( NewId ) {}

    /// Constructor using an array of nodes
    WaveEquationElement(IndexType NewId, const NodesArrayType& ThisNodes) : Element(NewId, ThisNodes) {}
    
    /// Constructor using Geometry
    WaveEquationElement(IndexType NewId, GeometryType::Pointer pGeometry) : Element( NewId, pGeometry ) {}
    
    /// Constructor using Properties
    WaveEquationElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : Element( NewId, pGeometry, pProperties ) 
    {
        mThisIntegrationMethod = this->GetGeometry().GetDefaultIntegrationMethod();
    }

    /// Destructor
    virtual ~WaveEquationElement() {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;
    
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const;
    
    int Check(const ProcessInfo& rCurrentProcessInfo);
              
    void GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo );
    
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   
	void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo );

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,ProcessInfo& rCurrentProcessInfo );
    
    void CalculateRightHandSide(VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo );
    
    void GetValuesVector(Vector& rValues, int Step = 0);
    
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0);
    
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0);
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
    
    /// Member Variables
    
    GeometryData::IntegrationMethod mThisIntegrationMethod;
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
	void CalculateAll(MatrixType& rLeftHandSideMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo );
	
	void CalculateLHS(MatrixType& rLeftHandSideMatrix,ProcessInfo& rCurrentProcessInfo );
	
	void CalculateRHS(VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo );
        
    void CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& detJ, const double& weight);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
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

    /// Assignment operator.
    WaveEquationElement & operator=(WaveEquationElement const& rOther);

    /// Copy constructor.
    WaveEquationElement(WaveEquationElement const& rOther);


}; // Class WaveEquationElement

} // namespace Kratos

#endif // KRATOS_WAVE_EQUATION_ELEMENT_H_INCLUDED  defined 
