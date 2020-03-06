//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#ifndef KRATOS_RANS_EVM_K_EPSILON_VMS_MONOLITHIC_WALL_H
#define KRATOS_RANS_EVM_K_EPSILON_VMS_MONOLITHIC_WALL_H

// System includes

// External includes

// Project includes

// Application includes
#include "custom_conditions/monolithic_wall_condition.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
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

/**
 * @brief
 *
 * @tparam TDim       Dimensionality of the condition (2D or 3D)
 * @tparam TNumNodes  Number of nodes in the condition
 */
template <unsigned int TDim, unsigned int TNumNodes = TDim>
class RansEvmKEpsilonVmsMonolithicWall : public MonolithicWallCondition<TDim, TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansEvmKEpsilonVmsMonolithicWall
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(RansEvmKEpsilonVmsMonolithicWall);

    using BaseType = MonolithicWallCondition<TDim, TNumNodes>;
    using NodeType = Node<3>;
    using PropertiesType = Properties;
    using GeometryType = Geometry<NodeType>;
    using NodesArrayType = Geometry<NodeType>::PointsArrayType;
    using VectorType = Vector;
    using MatrixType = Matrix;
    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new condition
      */
    explicit RansEvmKEpsilonVmsMonolithicWall(IndexType NewId = 0)
        : BaseType(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    RansEvmKEpsilonVmsMonolithicWall(IndexType NewId, const NodesArrayType& ThisNodes)
        : BaseType(NewId, ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    RansEvmKEpsilonVmsMonolithicWall(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    RansEvmKEpsilonVmsMonolithicWall(IndexType NewId,
                                     GeometryType::Pointer pGeometry,
                                     PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor.
    RansEvmKEpsilonVmsMonolithicWall(RansEvmKEpsilonVmsMonolithicWall const& rOther)
        : BaseType(rOther)
    {
    }

    /// Destructor.
    ~RansEvmKEpsilonVmsMonolithicWall() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    RansEvmKEpsilonVmsMonolithicWall& operator=(RansEvmKEpsilonVmsMonolithicWall const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /// Create a new RansEvmKEpsilonVmsMonolithicWall object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties) const override;

    Condition::Pointer Create(IndexType NewId,
                              GeometryType::Pointer pGeom,
                              PropertiesType::Pointer pProperties) const override;

    /**
     * Clones the selected element variables, creating a new one
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */

    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

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
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

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

    /**
     * @brief Applies wall law based on logarithmic wall law
     *
     * This applies wall law based on the following equation for the period of
     * time where RANS model is not activated for $y+ \leq 11.06$
     * \[
     *     u^+ = y^+
     * \]
     * for $y^+ > 11.06$
     * \[
     *     u^+ = \frac{1}{\kappa}ln\left(y^+\right\) + \beta
     * \]
     *
     * @param rLocalMatrix         Left hand side matrix
     * @param rLocalVector         Right hand side vector
     * @param rCurrentProcessInfo  Current process info from model part
     */
    virtual void ApplyLogarithmicWallLaw(MatrixType& rLocalMatrix,
                                         VectorType& rLocalVector,
                                         ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Applies rans based wall law
     *
     * This method calculate left hand side matrix and right hand side vector for following equation
     *
     * \[
     *      u_\tau = max\left(C_\mu^0.25 \sqrt{k}, \frac{||\underline{u}||}{\frac{1}{\kappa}ln(y^+)+\beta}\right)
     * \]
     *
     * integration point value = \rho \frac{u_\tau^2}{||\underline{u}||}\underline{u}
     *
     * @param rLocalMatrix         Left hand side matrix
     * @param rLocalVector         Right hand side vector
     * @param rCurrentProcessInfo  Current process info from model part
     */
    virtual void ApplyRansBasedWallLaw(MatrixType& rLocalMatrix,
                                       VectorType& rLocalVector,
                                       ProcessInfo& rCurrentProcessInfo);

    void ApplyWallLaw(MatrixType& rLocalMatrix,
                      VectorType& rLocalVector,
                      ProcessInfo& rCurrentProcessInfo) override;

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

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

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

}; // Class RansEvmKEpsilonVmsMonolithicWall

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::istream& operator>>(std::istream& rIStream,
                                RansEvmKEpsilonVmsMonolithicWall<TDim, TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansEvmKEpsilonVmsMonolithicWall<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_EVM_K_EPSILON_VMS_MONOLITHIC_WALL_H
