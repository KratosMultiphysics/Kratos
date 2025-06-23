// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ /
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/condition.h"
#include "geometries/coupling_geometry.h"

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

/**
 * @ingroup ContactStructuralMechanicsApplication
 * @class PairedCondition
 * @brief This is a base class for the conditions paired
 * @details This is a base class for the conditions paired, it is basically equal to the base condition, with a pointer to the paired geometry
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) PairedCondition
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of PairedCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( PairedCondition );

    /// Basetype definition
    using BaseType = Condition;

    /// Point type definition
    using PointType = Point;

    /// Geometry type definition
    using GeometryType = Geometry<Node>;

    /// Coupling geometry type definition
    using CouplingGeometryType = CouplingGeometry<Node>;

    /// Vector type definition
    using VectorType = BaseType::VectorType;

    /// Matrix type definition
    using MatrixType = BaseType::MatrixType;

    /// Index type definition
    using IndexType = BaseType::IndexType;

    /// Geometry pointer type definition
    using GeometryPointerType = BaseType::GeometryType::Pointer;

    /// Nodes array type definition
    using NodesArrayType = BaseType::NodesArrayType;

    /// Properties pointer type definition
    using PropertiesPointerType = BaseType::PropertiesType::Pointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    PairedCondition()
        : Condition()
    {}

    // Constructor 1
    PairedCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        ) :Condition(NewId, Kratos::make_shared<CouplingGeometryType>(pGeometry, nullptr))
    {
        KRATOS_WARNING_FIRST_N("PairedCondition", 10) << "This class pairs two geometries, please use the other constructor (the one with two geometries as input)" << std::endl;
    }

    // Constructor 2
    PairedCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        ) :Condition( NewId, Kratos::make_shared<CouplingGeometryType>(pGeometry, nullptr), pProperties )
    {
        KRATOS_WARNING_FIRST_N("PairedCondition", 10) << "This class pairs two geometries, please use the other constructor (the one with two geometries as input)" << std::endl;
    }

    // Constructor 3
    PairedCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties,
        GeometryType::Pointer pPairedGeometry
        )
        :Condition( NewId, Kratos::make_shared<CouplingGeometryType>(pGeometry, pPairedGeometry), pProperties )
    {}

    ///Copy constructor
    PairedCondition( PairedCondition const& rOther){}

    /// Destructor.
    ~PairedCondition() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Called at the beginning of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

   /**
    * @brief Called at the beginning of each solution step
    * @param rCurrentProcessInfo the current process info instance
    */
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

   /**
    * @brief Called at the beginning of each iteration
    * @param rCurrentProcessInfo the current process info instance
    */
    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Creates a new element pointer from an array of nodes
     * @param NewId the ID of the new element
     * @param rThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& rThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new element pointer from an existing geometry
     * @param NewId the ID of the new element
     * @param pGeometry the  geometry taken to create the condition
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new element pointer from an existing geometry
     * @param NewId the ID of the new element
     * @param pGeometry the  geometry taken to create the condition
     * @param pProperties the properties assigned to the new element
     * @param pPairedGeom the paired geometry
     * @return a Pointer to the new element
     */
    virtual Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties,
        GeometryType::Pointer pPairedGeom
        ) const;

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This method returns the parent geometry
     * @return The slave geometry (slave in the definition of Popp which is the opposite of the standard)
     */
    GeometryType::Pointer pGetParentGeometry()
    {
        return this->GetGeometry().pGetGeometryPart(CouplingGeometryType::Master);
    }

    /**
     * @brief This method returns the parent geometry (constant version)
     * @return The slave geometry (slave in the definition of Popp which is the opposite of the standard)
     */
    GeometryType::Pointer const pGetParentGeometry() const
    {
        return this->GetGeometry().pGetGeometryPart(CouplingGeometryType::Master);
    }

    /**
     * @brief This method returns the paired geometry
     * @return The master geometry (master in the definition of Popp which is the opposite of the standard)
     */
    GeometryType::Pointer pGetPairedGeometry()
    {
        return this->GetGeometry().pGetGeometryPart(CouplingGeometryType::Slave);
    }

    /**
     * @brief This method returns the paired geometry (constant version)
     * @return The master geometry (master in the definition of Popp which is the opposite of the standard)
     */
    GeometryType::Pointer const pGetPairedGeometry() const
    {
        return this->GetGeometry().pGetGeometryPart(CouplingGeometryType::Slave);
    }

    /**
     * @brief This method returns the parent geometry
     * @return The slave geometry (slave in the definition of Popp which is the opposite of the standard)
     */
    GeometryType& GetParentGeometry()
    {
        return this->GetGeometry().GetGeometryPart(CouplingGeometryType::Master);
    }

    /**
     * @brief This method returns the parent geometry (constant version)
     * @return The slave geometry (slave in the definition of Popp which is the opposite of the standard)
     */
    GeometryType const& GetParentGeometry() const
    {
        return this->GetGeometry().GetGeometryPart(CouplingGeometryType::Master);
    }

    /**
     * @brief This method returns the paired geometry
     * @return The master geometry (master in the definition of Popp which is the opposite of the standard)
     */
    GeometryType& GetPairedGeometry()
    {
        return this->GetGeometry().GetGeometryPart(CouplingGeometryType::Slave);
    }

    /**
     * @brief This method returns the paired geometry (constant version)
     * @return The master geometry (master in the definition of Popp which is the opposite of the standard)
     */
    GeometryType const& GetPairedGeometry() const
    {
        return this->GetGeometry().GetGeometryPart(CouplingGeometryType::Slave);
    }

    /**
     * @brief This method sets the paired normal
     * @param rPairedNormal The master geometry normal
     */
    void SetPairedNormal(const array_1d<double, 3>& rPairedNormal)
    {
        noalias(mPairedNormal) = rPairedNormal;
    }

    /**
     * @brief This method returns the paired normal
     * @return The master geometry normal
     */
    array_1d<double, 3> const& GetPairedNormal() const
    {
        return mPairedNormal;
    }

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
        buffer << "PairedCondition #" << this->Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PairedCondition #" << this->Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        PrintInfo(rOStream);
        this->GetParentGeometry().PrintData(rOStream);
        this->GetPairedGeometry().PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    array_1d<double, 3> mPairedNormal = ZeroVector(3);

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

    // Serialization

    friend class Serializer;

    /**
     * @brief Saves the PairedCondition object to a serializer.
     * @param rSerializer the serializer to save to
     */
    void save(Serializer& rSerializer) const override;

    /**
     * @brief Loads the PairedCondition from a serializer.
     * @param rSerializer the serializer to load from
     */
    void load(Serializer& rSerializer) override;

    ///@}

}; // Class PairedCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.
