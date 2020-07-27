//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Mayu Sakuma
//

#if !defined(KRATOS_KINEMATIC_SIMULATION_SPECTRAL_CONSTANT_V_H_INCLUDED)
#define KRATOS_KINEMATIC_SIMULATION_SPECTRAL_CONSTANT_V_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/kinematic_simulation/kinematic_simulation_element.h"
#include "includes/checks.h"
#include "includes/define.h"

// Application includes
#include "includes/variables.h"
#include "rans_application_variables.h"

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

template <unsigned int TNumNodes>
class KinematicSimulationSpectralConstantVElement : public KinematicSimulationElement<TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = KinematicSimulationElement<TNumNodes>;

    /// Node type (default is: Node<3>)
    using NodeType = Node<3>;

    /// Geometry type (using with given NodeType)
    using GeometryType = Geometry<NodeType>;

    using PropertiesType = Properties;

    /// Definition of nodes container type, redefined from GeometryType
    using NodesArrayType = GeometryType::PointsArrayType;

    /// Vector type for local contributions to the linear system
    using VectorType = Element::VectorType;

    /// Matrix type for local contributions to the linear system
    using MatrixType = Element::MatrixType;

    using IndexType = std::size_t;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of LaplaceElement
    KRATOS_CLASS_POINTER_DEFINITION(KinematicSimulationSpectralConstantVElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit KinematicSimulationSpectralConstantVElement(IndexType NewId = 0)
        : BaseType(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    KinematicSimulationSpectralConstantVElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : BaseType(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    KinematicSimulationSpectralConstantVElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    KinematicSimulationSpectralConstantVElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    KinematicSimulationSpectralConstantVElement(KinematicSimulationSpectralConstantVElement const& rOther)
        : BaseType(rOther)
    {
    }

    /**
     * Destructor
     */
    ~KinematicSimulationSpectralConstantVElement() override = default;

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
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY;
        return Kratos::make_intrusive<KinematicSimulationSpectralConstantVElement<TNumNodes>>(
            NewId, Element::GetGeometry().Create(ThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY;
        return Kratos::make_intrusive<KinematicSimulationSpectralConstantVElement<TNumNodes>>(
            NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
    {
        KRATOS_TRY;
        return Kratos::make_intrusive<KinematicSimulationSpectralConstantVElement<TNumNodes>>(
            NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
        KRATOS_CATCH("");
    }

    const Variable<double>& GetVariable() const override
    {
        return SPECTRAL_CONSTANT_V;
    }

    /**
     * ELEMENTS inherited from this class have to implement next
     * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide
     * methods they can be managed internally with a private method to do the
     * same calculations only once: MANDATORY
     */


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
        buffer << "KinematicSimulationSpectralConstantVElement #" << this->Id();
        return buffer.str();
    }
    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "KinematicSimulationSpectralConstantVElement #" << this->Id();
    }

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
        KRATOS_TRY

        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);

        KRATOS_CATCH("");
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_TRY

        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);

        KRATOS_CATCH("");
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

}; // Class LaplaceElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

template <unsigned int TNumNodes>
inline std::istream& operator>>(std::istream& rIStream,
                                KinematicSimulationSpectralConstantVElement<TNumNodes>& rThis);

/// output stream function
template <unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const KinematicSimulationSpectralConstantVElement<TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_LAPLACE_ELEMENT_H_INCLUDED defined
