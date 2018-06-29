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

#if !defined(KRATOS_PRIMITIVE_VAR_ELEM_H_INCLUDED)
#define  KRATOS_PRIMITIVE_VAR_ELEM_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"

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

/// Implementation of a linear element for shallow water problems
template< unsigned int TNumNodes >
class PrimitiveVarElement : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of PrimitiveVarElement
    KRATOS_CLASS_POINTER_DEFINITION( PrimitiveVarElement );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PrimitiveVarElement() :
        Element()
    {}

    /// Constructor using a Geometry instance
    PrimitiveVarElement(IndexType NewId, GeometryType::Pointer pGeometry) :
        Element(NewId, pGeometry)
    {}

    /// Constructor using geometry and properties
    PrimitiveVarElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
        Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~ PrimitiveVarElement() {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Create a new Primitive variables element and return a pointer to it
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Element::Pointer(new PrimitiveVarElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("")
    }


    /// Check that all required data containers are properly initialized and registered in Kratos
    /**
     * @return 0 if no errors are detected.
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    /// Fill given vector with the linear system row index for the element's degrees of freedom
    /**
     * @param rResult
     * @param rCurrentProcessInfo
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    /// Fill given array with containing the element's degrees of freedom
    /**
     * @param rElementalDofList
     * @param rCurrentProcessInfo
     */
    void GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo) override;

    /// Evaluate the elemental contribution to the problem for turbulent viscosity.
    /**
     * @param rLeftHandSideMatrix Elemental left hand side matrix
     * @param rRightHandSideVector Elemental right hand side vector
     * @param rCurrentProcessInfo Reference to the ProcessInfo from the ModelPart containg the element
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

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

    struct ElementVariables
    {
        double dt_inv;
        double lumping_factor;
        double dyn_tau;
        double gravity;
        double manning2;
        double height_units;

        double height;
        array_1d<double,2> velocity;
        array_1d<double,2> momentum;
        array_1d<double,2> height_grad;
        BoundedMatrix<double,2,2> velocity_grad;
        double velocity_div;

        array_1d<double, TNumNodes*3> depth;
        array_1d<double, TNumNodes*3> rain;
        array_1d<double, TNumNodes*3> unknown;
        array_1d<double, TNumNodes*3> proj_unk;
    };

    void InitializeElement(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    void CalculateGeometry(BoundedMatrix<double, TNumNodes, 2>& rDN_DX, double& rArea);

    double ComputeElemSize(const BoundedMatrix<double, TNumNodes, 2>& rDN_DX);

    void GetNodalValues(ElementVariables& rVariables);

    void GetElementValues(const BoundedMatrix<double,TNumNodes, 2>& rDN_DX, ElementVariables& rVariables);

    void ComputeStabilizationParameters(const ElementVariables& rVariables,
                                        const double& rElemSize,
                                        double& rTauU,
                                        double& rTauH,
                                        double& rKdc);

    void ComputeAuxMatrices(
            const BoundedMatrix<double,TNumNodes, TNumNodes>& rNcontainer,
            const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
            const ElementVariables& rVariables,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rMassMatrixScalar,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rMassMatrixVector,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rScalarGrad,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rVectorDiv,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rScalarDiff,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rVectorDiff );

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

}; // Class PrimitiveVarElement

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_PRIMITIVE_VAR_ELEM_H_INCLUDED  defined
