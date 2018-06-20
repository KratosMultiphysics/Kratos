//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//

#if !defined(KRATOS_SHALLOW_ELEMENT_H_INCLUDED )
#define KRATOS_SHALLOW_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"

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

class ShallowElement : public Element
{
public:

    ///@name Type Definitions
    ///@{

    struct ElementData
    {
        double dt_inv;
        double lumping_factor;
        double c_tau;
        double gravity;
        double manning2;

        array_1d<double, 9> depth;
        array_1d<double, 9> rain;
        array_1d<double, 9> unknown;
        array_1d<double, 9> proj_unk;
    };

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of ShallowElement
    KRATOS_CLASS_POINTER_DEFINITION(ShallowElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    ShallowElement(IndexType NewId = 0)
    : Element(NewId)
    {}

    /**
     * Constructor using an array of nodes
     */
    ShallowElement(IndexType NewId, const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes)
    {}

    /**
     * Constructor using Geometry
     */
    ShallowElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}

    /**
     * Constructor using Properties
     */
    ShallowElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /**
     * Destructor
     */
    ~ShallowElement(){};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * ELEMENTS inherited from this class have to implement next
     * Create and Clone methods: MANDATORY
     */

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared<ShallowElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared<ShallowElement>(NewId, pGeom, pProperties);
    }

    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult: the elemental equation ID vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) override;

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList: the list of DOFs
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix only
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override;

    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the constitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Access for variables on Integration points.
     * This gives access to variables stored in the constitutive law on each integration point.
     * Specializations of element must specify the actual interface to the integration points!
     * Note, that these functions expect a std::vector of values for the specified variable type that
     * contains a value for each integration point!
     * GetValueOnIntegrationPoints: get the values for given Variable.
     * @param rVariable: the specified variable
     * @param rValues: where to store the values for the specified variable type at each integration point
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

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
        return "Shallow water element";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info() << Id();
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{

    static constexpr unsigned int msNodes = 3;
    static constexpr unsigned int msElemSize = 3 * msNodes;

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void InitializeElement(
        ElementData& rData,
        const ProcessInfo& rCurrentProcessInfo);

    void ComputeMassMatrices(
        const ElementData& rData,
        BoundedMatrix<double,msElemSize,msElemSize>& rVelMatrix,
        BoundedMatrix<double,msElemSize,msElemSize>& rHeightMatrix);

    void ComputeElementValues(
        const ElementData& rData,
        array_1d<double,2>& rVel,
        double& rHeight);

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
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

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

}; // Class ShallowElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_SHALLOW_ELEMENT_H_INCLUDED  defined
