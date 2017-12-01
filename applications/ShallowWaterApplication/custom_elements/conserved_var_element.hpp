//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#if !defined(KRATOS_CONSERVED_VAR_ELEM_H_INCLUDED)
#define  KRATOS_CONSERVED_VAR_ELEM_H_INCLUDED 

// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/serializer.h"

namespace Kratos
{
///@addtogroup ShallowWaterApplication
///@{

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

/// Implementation of a linear element for shallow water interpolating coserved variables
template< unsigned int TNumNodes >
class ConservedVarElement : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of ConservedVarElement
    KRATOS_CLASS_POINTER_DEFINITION( ConservedVarElement );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ConservedVarElement() :
        Element()
    {}

    /// Constructor using a Geometry instance
    ConservedVarElement(IndexType NewId, GeometryType::Pointer pGeometry) :
        Element(NewId, pGeometry)
    {}

    /// Constructor using geometry and properties
    ConservedVarElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
        Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~ ConservedVarElement() {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Create a new Primitive variables element and return a pointer to it
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
        return Element::Pointer(new ConservedVarElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("")
    }

    /// Check that all required data containers are properly initialized and registered in Kratos
    /** 
     * @return 0 if no errors are detected.
     */
    int Check(const ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    void GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo);

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


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

    void CalculateGeometry(boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX, double& rArea);

    double ComputeElemSize(boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX);

    void GetNodalValues(array_1d<double,TNumNodes*3>& rDepth, array_1d<double,TNumNodes*3>& rRain, array_1d<double,TNumNodes*3>& rUnkn, array_1d<double,TNumNodes*3>& rProj);

    void GetElementValues(const boost::numeric::ublas::bounded_matrix<double,TNumNodes,2>& rDN_DX, const array_1d<double,TNumNodes*3>& rNodalVar, array_1d<double,2>& rMomentum, double& rDivU, double& rHeight, array_1d<double,2>& rHeightGrad);

    void ComputeStabilizationParameters(const double& rHeight, const array_1d<double,2>& rHeightGrad, const double& rElemSize, double& rTauU, double& rTauH, double& rKdc, const ProcessInfo& rCurrentProcessInfo);

    void CalculateLumpedMassMatrix(boost::numeric::ublas::bounded_matrix<double, TNumNodes*3, TNumNodes*3>& rMassMatrix);

    double mGravity;
    double mHeightUnitConvert;

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
    ///@name Serialization
    ///@{

    friend class Serializer;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

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


}; // Class ConservedVarElement

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CONSERVED_VAR_ELEM_H_INCLUDED  defined
