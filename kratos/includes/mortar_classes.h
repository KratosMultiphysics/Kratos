//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "utilities/mortar_utilities.h"
#include "includes/variables.h"
#include "includes/mapping_variables.h"
#include "includes/serializer.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The definition of the size type
    using SizeType = std::size_t;

///@}
///@name  Enum's
///@{
    /**
    * @brief This enum defines a "hash" used to identify in which combination of cuts the point is found when the mortar segmentation is performed
    */
    enum class PointBelongs
    {
        Master = 0,      /// The point belongs to the master entity
        Slave = 1,       /// The point belongs to the slave entity
        Intersection = 2 /// The point belongs to an intersection of entities
    };

    /**
    * @brief Case for 2D line intersected with another 2D line
    */
    enum class PointBelongsLine2D2N
    {
        SlaveLine2D2N0 = 0,       /// The point belongs to the first slave node of a 2D line element
        SlaveLine2D2N1 = 1,       /// The point belongs to the second slave node of a 2D line element
        MasterLine2D2N0 = 2,      /// The point belongs to the first master node of a 2D line element
        MasterLine2D2N1 = 3,      /// The point belongs to the second master node of a 2D line element
        IntersectionLine2D2N = 4  /// The point belongs to an intersection of nodes in a 2D line element
    };

    /**
    * @brief Case for 3D triangle intersected with another 3D triangle
    */
    enum class PointBelongsTriangle3D3N
    {
        SlaveTriangle3D3N0           = 0,
        SlaveTriangle3D3N1           = 1,
        SlaveTriangle3D3N2           = 2,
        MasterTriangle3D3N0          = 3,
        MasterTriangle3D3N1          = 4,
        MasterTriangle3D3N2          = 5,
        IntersectionTriangle3D3N     = 6,
        IntersectionTriangle3D3N0101 = 10106,
        IntersectionTriangle3D3N1001 = 10016,
        IntersectionTriangle3D3N1201 = 10216,
        IntersectionTriangle3D3N2101 = 10126,
        IntersectionTriangle3D3N0201 = 10206,
        IntersectionTriangle3D3N2001 = 10026,
        IntersectionTriangle3D3N0110 = 1106,
        IntersectionTriangle3D3N1010 = 1016,
        IntersectionTriangle3D3N1210 = 1216,
        IntersectionTriangle3D3N2110 = 1126,
        IntersectionTriangle3D3N0210 = 1206,
        IntersectionTriangle3D3N2010 = 1026,
        IntersectionTriangle3D3N0112 = 21106,
        IntersectionTriangle3D3N1012 = 21016,
        IntersectionTriangle3D3N1212 = 21216,
        IntersectionTriangle3D3N2112 = 21126,
        IntersectionTriangle3D3N0212 = 21206,
        IntersectionTriangle3D3N2012 = 21026,
        IntersectionTriangle3D3N0121 = 12106,
        IntersectionTriangle3D3N1021 = 12016,
        IntersectionTriangle3D3N1221 = 12216,
        IntersectionTriangle3D3N2121 = 12126,
        IntersectionTriangle3D3N0221 = 12206,
        IntersectionTriangle3D3N2021 = 12026,
        IntersectionTriangle3D3N0102 = 20106,
        IntersectionTriangle3D3N1002 = 20016,
        IntersectionTriangle3D3N1202 = 20216,
        IntersectionTriangle3D3N2102 = 20126,
        IntersectionTriangle3D3N0202 = 20206,
        IntersectionTriangle3D3N2002 = 20026,
        IntersectionTriangle3D3N0120 = 2106,
        IntersectionTriangle3D3N1020 = 2016,
        IntersectionTriangle3D3N1220 = 2216,
        IntersectionTriangle3D3N2120 = 2126,
        IntersectionTriangle3D3N0220 = 2206,
        IntersectionTriangle3D3N2020 = 2026
    };

    /**
    * @brief Case for 3D quadrilateral intersected with another 3D quadrilateral
    */
    enum class PointBelongsQuadrilateral3D4N
    {
        SlaveQuadrilateral3D4N0           = 0,
        SlaveQuadrilateral3D4N1           = 1,
        SlaveQuadrilateral3D4N2           = 2,
        SlaveQuadrilateral3D4N3           = 3,
        MasterQuadrilateral3D4N0          = 4,
        MasterQuadrilateral3D4N1          = 5,
        MasterQuadrilateral3D4N2          = 6,
        MasterQuadrilateral3D4N3          = 7,
        IntersectionQuadrilateral3D4N     = 8,
        IntersectionQuadrilateral3D4N0101 = 10108,
        IntersectionQuadrilateral3D4N1001 = 10018,
        IntersectionQuadrilateral3D4N1201 = 10218,
        IntersectionQuadrilateral3D4N2101 = 10128,
        IntersectionQuadrilateral3D4N2301 = 10328,
        IntersectionQuadrilateral3D4N3201 = 10238,
        IntersectionQuadrilateral3D4N3001 = 10038,
        IntersectionQuadrilateral3D4N0301 = 10308,
        IntersectionQuadrilateral3D4N0110 = 1108,
        IntersectionQuadrilateral3D4N1010 = 1018,
        IntersectionQuadrilateral3D4N1210 = 1218,
        IntersectionQuadrilateral3D4N2110 = 1128,
        IntersectionQuadrilateral3D4N2310 = 1328,
        IntersectionQuadrilateral3D4N3210 = 1238,
        IntersectionQuadrilateral3D4N3010 = 1038,
        IntersectionQuadrilateral3D4N0310 = 1308,
        IntersectionQuadrilateral3D4N0112 = 21108,
        IntersectionQuadrilateral3D4N1012 = 21018,
        IntersectionQuadrilateral3D4N1212 = 21218,
        IntersectionQuadrilateral3D4N2112 = 21128,
        IntersectionQuadrilateral3D4N2312 = 21328,
        IntersectionQuadrilateral3D4N3212 = 21238,
        IntersectionQuadrilateral3D4N3012 = 21038,
        IntersectionQuadrilateral3D4N0312 = 21308,
        IntersectionQuadrilateral3D4N0121 = 12108,
        IntersectionQuadrilateral3D4N1021 = 12018,
        IntersectionQuadrilateral3D4N1221 = 12218,
        IntersectionQuadrilateral3D4N2121 = 12128,
        IntersectionQuadrilateral3D4N2321 = 12328,
        IntersectionQuadrilateral3D4N3221 = 12238,
        IntersectionQuadrilateral3D4N3021 = 12038,
        IntersectionQuadrilateral3D4N0321 = 12308,
        IntersectionQuadrilateral3D4N0123 = 32108,
        IntersectionQuadrilateral3D4N1023 = 32018,
        IntersectionQuadrilateral3D4N1223 = 32218,
        IntersectionQuadrilateral3D4N2123 = 32128,
        IntersectionQuadrilateral3D4N2323 = 32328,
        IntersectionQuadrilateral3D4N3223 = 32238,
        IntersectionQuadrilateral3D4N3023 = 32038,
        IntersectionQuadrilateral3D4N0323 = 32308,
        IntersectionQuadrilateral3D4N0132 = 23108,
        IntersectionQuadrilateral3D4N1032 = 23018,
        IntersectionQuadrilateral3D4N1232 = 23218,
        IntersectionQuadrilateral3D4N2132 = 23128,
        IntersectionQuadrilateral3D4N2332 = 23328,
        IntersectionQuadrilateral3D4N3232 = 23238,
        IntersectionQuadrilateral3D4N3032 = 23038,
        IntersectionQuadrilateral3D4N0332 = 23308,
        IntersectionQuadrilateral3D4N0130 = 3108,
        IntersectionQuadrilateral3D4N1030 = 3018,
        IntersectionQuadrilateral3D4N1230 = 3218,
        IntersectionQuadrilateral3D4N2130 = 3128,
        IntersectionQuadrilateral3D4N2330 = 3328,
        IntersectionQuadrilateral3D4N3230 = 3238,
        IntersectionQuadrilateral3D4N3030 = 3038,
        IntersectionQuadrilateral3D4N0330 = 3308,
        IntersectionQuadrilateral3D4N0103 = 30108,
        IntersectionQuadrilateral3D4N1003 = 30018,
        IntersectionQuadrilateral3D4N1203 = 30218,
        IntersectionQuadrilateral3D4N2103 = 30128,
        IntersectionQuadrilateral3D4N2303 = 30328,
        IntersectionQuadrilateral3D4N3203 = 30238,
        IntersectionQuadrilateral3D4N3003 = 30038,
        IntersectionQuadrilateral3D4N0303 = 30308
    };

    /**
    * @brief Case for 3D triangle intersected with quadrilateral 3D
    */
    enum class PointBelongsQuadrilateral3D4NTriangle3D3N
    {
        SlaveQuadrilateral3D4N0                       = 0,
        SlaveQuadrilateral3D4N1                       = 1,
        SlaveQuadrilateral3D4N2                       = 2,
        SlaveQuadrilateral3D4N3                       = 3,
        MasterTriangle3D3N0                           = 4,
        MasterTriangle3D3N1                           = 5,
        MasterTriangle3D3N2                           = 6,
        IntersectionQuadrilateral3D4NTriangle3D3N     = 7,
        IntersectionQuadrilateral3D4NTriangle3D3N0101 = 10107,
        IntersectionQuadrilateral3D4NTriangle3D3N1001 = 1107,
        IntersectionQuadrilateral3D4NTriangle3D3N1201 = 21107,
        IntersectionQuadrilateral3D4NTriangle3D3N2101 = 12107,
        IntersectionQuadrilateral3D4NTriangle3D3N2301 = 32107,
        IntersectionQuadrilateral3D4NTriangle3D3N3201 = 23107,
        IntersectionQuadrilateral3D4NTriangle3D3N3001 = 3107,
        IntersectionQuadrilateral3D4NTriangle3D3N0301 = 30107,
        IntersectionQuadrilateral3D4NTriangle3D3N0110 = 10017,
        IntersectionQuadrilateral3D4NTriangle3D3N1010 = 1017,
        IntersectionQuadrilateral3D4NTriangle3D3N1210 = 21017,
        IntersectionQuadrilateral3D4NTriangle3D3N2110 = 12017,
        IntersectionQuadrilateral3D4NTriangle3D3N2310 = 32017,
        IntersectionQuadrilateral3D4NTriangle3D3N3210 = 23017,
        IntersectionQuadrilateral3D4NTriangle3D3N3010 = 3017,
        IntersectionQuadrilateral3D4NTriangle3D3N0310 = 30017,
        IntersectionQuadrilateral3D4NTriangle3D3N0112 = 10217,
        IntersectionQuadrilateral3D4NTriangle3D3N1012 = 1217,
        IntersectionQuadrilateral3D4NTriangle3D3N1212 = 21217,
        IntersectionQuadrilateral3D4NTriangle3D3N2112 = 12217,
        IntersectionQuadrilateral3D4NTriangle3D3N2312 = 32217,
        IntersectionQuadrilateral3D4NTriangle3D3N3212 = 23217,
        IntersectionQuadrilateral3D4NTriangle3D3N3012 = 3217,
        IntersectionQuadrilateral3D4NTriangle3D3N0312 = 30217,
        IntersectionQuadrilateral3D4NTriangle3D3N0121 = 10127,
        IntersectionQuadrilateral3D4NTriangle3D3N1021 = 1127,
        IntersectionQuadrilateral3D4NTriangle3D3N1221 = 21127,
        IntersectionQuadrilateral3D4NTriangle3D3N2121 = 12127,
        IntersectionQuadrilateral3D4NTriangle3D3N2321 = 32127,
        IntersectionQuadrilateral3D4NTriangle3D3N3221 = 23127,
        IntersectionQuadrilateral3D4NTriangle3D3N3021 = 3127,
        IntersectionQuadrilateral3D4NTriangle3D3N0321 = 30127,
        IntersectionQuadrilateral3D4NTriangle3D3N0120 = 10027,
        IntersectionQuadrilateral3D4NTriangle3D3N1020 = 1027,
        IntersectionQuadrilateral3D4NTriangle3D3N1220 = 21027,
        IntersectionQuadrilateral3D4NTriangle3D3N2120 = 12027,
        IntersectionQuadrilateral3D4NTriangle3D3N2320 = 32027,
        IntersectionQuadrilateral3D4NTriangle3D3N3220 = 23027,
        IntersectionQuadrilateral3D4NTriangle3D3N3020 = 3027,
        IntersectionQuadrilateral3D4NTriangle3D3N0320 = 30027,
        IntersectionQuadrilateral3D4NTriangle3D3N0102 = 10207,
        IntersectionQuadrilateral3D4NTriangle3D3N1002 = 1207,
        IntersectionQuadrilateral3D4NTriangle3D3N1202 = 21207,
        IntersectionQuadrilateral3D4NTriangle3D3N2102 = 12207,
        IntersectionQuadrilateral3D4NTriangle3D3N2302 = 32207,
        IntersectionQuadrilateral3D4NTriangle3D3N3202 = 23207,
        IntersectionQuadrilateral3D4NTriangle3D3N3002 = 3207,
        IntersectionQuadrilateral3D4NTriangle3D3N0302 = 30207
    };

    /**
    * @brief Case for 3D quadrilateral intersected with triangle 3D
    */
    enum class PointBelongsTriangle3D3NQuadrilateral3D4N
    {
        SlaveTriangle3D3N0                            = 0,
        SlaveTriangle3D3N1                            = 1,
        SlaveTriangle3D3N2                            = 2,
        MasterQuadrilateral3D4N0                      = 3,
        MasterQuadrilateral3D4N1                      = 4,
        MasterQuadrilateral3D4N2                      = 5,
        MasterQuadrilateral3D4N3                      = 6,
        IntersectionTriangle3D3NQuadrilateral3D4N     = 7,
        IntersectionTriangle3D3NQuadrilateral3D4N0101 = 10107,
        IntersectionTriangle3D3NQuadrilateral3D4N0110 = 1107,
        IntersectionTriangle3D3NQuadrilateral3D4N0112 = 21107,
        IntersectionTriangle3D3NQuadrilateral3D4N0121 = 12107,
        IntersectionTriangle3D3NQuadrilateral3D4N0123 = 32107,
        IntersectionTriangle3D3NQuadrilateral3D4N0132 = 23107,
        IntersectionTriangle3D3NQuadrilateral3D4N0130 = 3107,
        IntersectionTriangle3D3NQuadrilateral3D4N0103 = 30107,
        IntersectionTriangle3D3NQuadrilateral3D4N1001 = 10017,
        IntersectionTriangle3D3NQuadrilateral3D4N1010 = 1017,
        IntersectionTriangle3D3NQuadrilateral3D4N1012 = 21017,
        IntersectionTriangle3D3NQuadrilateral3D4N1021 = 12017,
        IntersectionTriangle3D3NQuadrilateral3D4N1023 = 32017,
        IntersectionTriangle3D3NQuadrilateral3D4N1032 = 23017,
        IntersectionTriangle3D3NQuadrilateral3D4N1030 = 3017,
        IntersectionTriangle3D3NQuadrilateral3D4N1003 = 30017,
        IntersectionTriangle3D3NQuadrilateral3D4N1201 = 10217,
        IntersectionTriangle3D3NQuadrilateral3D4N1210 = 1217,
        IntersectionTriangle3D3NQuadrilateral3D4N1212 = 21217,
        IntersectionTriangle3D3NQuadrilateral3D4N1221 = 12217,
        IntersectionTriangle3D3NQuadrilateral3D4N1223 = 32217,
        IntersectionTriangle3D3NQuadrilateral3D4N1232 = 23217,
        IntersectionTriangle3D3NQuadrilateral3D4N1230 = 3217,
        IntersectionTriangle3D3NQuadrilateral3D4N1203 = 30217,
        IntersectionTriangle3D3NQuadrilateral3D4N2101 = 10127,
        IntersectionTriangle3D3NQuadrilateral3D4N2110 = 1127,
        IntersectionTriangle3D3NQuadrilateral3D4N2112 = 21127,
        IntersectionTriangle3D3NQuadrilateral3D4N2121 = 12127,
        IntersectionTriangle3D3NQuadrilateral3D4N2123 = 32127,
        IntersectionTriangle3D3NQuadrilateral3D4N2132 = 23127,
        IntersectionTriangle3D3NQuadrilateral3D4N2130 = 3127,
        IntersectionTriangle3D3NQuadrilateral3D4N2103 = 30127,
        IntersectionTriangle3D3NQuadrilateral3D4N2001 = 10027,
        IntersectionTriangle3D3NQuadrilateral3D4N2010 = 1027,
        IntersectionTriangle3D3NQuadrilateral3D4N2012 = 21027,
        IntersectionTriangle3D3NQuadrilateral3D4N2021 = 12027,
        IntersectionTriangle3D3NQuadrilateral3D4N2023 = 32027,
        IntersectionTriangle3D3NQuadrilateral3D4N2032 = 23027,
        IntersectionTriangle3D3NQuadrilateral3D4N2030 = 3027,
        IntersectionTriangle3D3NQuadrilateral3D4N2003 = 30027,
        IntersectionTriangle3D3NQuadrilateral3D4N0201 = 10207,
        IntersectionTriangle3D3NQuadrilateral3D4N0210 = 1207,
        IntersectionTriangle3D3NQuadrilateral3D4N0212 = 21207,
        IntersectionTriangle3D3NQuadrilateral3D4N0221 = 12207,
        IntersectionTriangle3D3NQuadrilateral3D4N0223 = 32207,
        IntersectionTriangle3D3NQuadrilateral3D4N0232 = 23207,
        IntersectionTriangle3D3NQuadrilateral3D4N0230 = 3207,
        IntersectionTriangle3D3NQuadrilateral3D4N0203 = 30207
    };

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class MortarKinematicVariables
 * @ingroup KratosCore
 * @brief MortarKinematicVariables
 * @details This is the definition of the kinematic variables used on the mortar operators assemble, which means three shape functions (one for the slave , one for the master and the third for the Lagrange Multipliers), and the jacobian in the corresponding Gauss point
 * @author Vicente Mataix Ferrandiz
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template< const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
class MortarKinematicVariables
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MortarKinematicVariables
    KRATOS_CLASS_POINTER_DEFINITION( MortarKinematicVariables );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    MortarKinematicVariables() = default;

    /// Destructor
    virtual ~MortarKinematicVariables() = default;

    ///@}
    ///@name Member Variables
    ///@{

    // Shape functions for contact pair
    Vector NMaster = Vector(TNumNodesMaster, 0.0);
    Vector NSlave = Vector(TNumNodes, 0.0);
    Vector PhiLagrangeMultipliers = Vector(TNumNodes, 0.0);

    // Determinant of slave cell's jacobian
    double DetjSlave = 0.0;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * This method initialized the operators
     */
    void Initialize()
    {
        // Shape functions
        noalias(NMaster)                = ZeroVector(TNumNodesMaster);
        noalias(NSlave)                 = ZeroVector(TNumNodes);
        noalias(PhiLagrangeMultipliers) = ZeroVector(TNumNodes);

        // Jacobian of slave
        DetjSlave = 0.0;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "MortarKinematicVariables with " << TNumNodes << " nodes (slave) and " << TNumNodesMaster << " nodes (master)";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MortarKinematicVariables with " << TNumNodes << " nodes (slave) and " << TNumNodesMaster << " nodes (master)" << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        // Print the base class data
        rOStream << "NSlave:" << NSlave << "\n";
        rOStream << "NMaster: " <<NMaster << "\n";
        rOStream << "PhiLagrangeMultipliers: "<< PhiLagrangeMultipliers << "\n";
        rOStream << "DetjSlave: " << DetjSlave << std::endl;
    }

    ///@}
private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        rSerializer.save("NMaster", NMaster);
        rSerializer.save("NSlave", NSlave);
        rSerializer.save("PhiLagrangeMultipliers", PhiLagrangeMultipliers);
        rSerializer.save("DetjSlave", DetjSlave);
    }

    virtual void load(Serializer& rSerializer)
    {
        rSerializer.load("NMaster", NMaster);
        rSerializer.load("NSlave", NSlave);
        rSerializer.load("PhiLagrangeMultipliers", PhiLagrangeMultipliers);
        rSerializer.load("DetjSlave", DetjSlave);
    }

    ///@}
}; // Class MortarKinematicVariables

/**
 * @class MortarKinematicVariablesWithDerivatives
 * @ingroup KratosCore
 * @brief MortarKinematicVariablesWithDerivatives
 * @details This class derives from MortarKinematicVariables and includes additionally to the variables of the previous class, the variables needed to define the directional derivatives of the mortar operators, like the gradients of the shape functions and the jacobians
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template< const SizeType TDim, const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
class MortarKinematicVariablesWithDerivatives
    : public MortarKinematicVariables<TNumNodes, TNumNodesMaster>
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the base class
    using BaseClassType = MortarKinematicVariables<TNumNodes, TNumNodesMaster>;

    /// Counted pointer of MortarKinematicVariables
    KRATOS_CLASS_POINTER_DEFINITION( MortarKinematicVariablesWithDerivatives );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    MortarKinematicVariablesWithDerivatives()= default;

    /// Destructor
    ~MortarKinematicVariablesWithDerivatives() override= default;

    ///@}
    ///@name Member Variables
    ///@{

    // Shape functions local derivatives for contact pair
    Matrix DNDeMaster = ScalarMatrix(TNumNodesMaster, TDim - 1, 0.0);
    Matrix DNDeSlave = ScalarMatrix(TNumNodes, TDim - 1, 0.0);

    /*
    * Jacobians in current configuration on all integration points of slave segment
    * Only those two variables contain info on all GP
    * other variables contain info only on the currently-calculated GP
    */
    Matrix jSlave = ScalarMatrix(TDim, TDim - 1, 0.0);
    Matrix jMaster = ScalarMatrix(TDim, TDim - 1, 0.0);

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method initialized the operators
     */
    void Initialize()
    {
        BaseClassType::Initialize();

        // Shape functions local derivatives
        noalias(DNDeMaster) = ZeroMatrix(TNumNodesMaster, TDim - 1);
        noalias(DNDeSlave)  = ZeroMatrix(TNumNodes, TDim - 1);

        // Jacobians on all integration points
        noalias(jSlave)  = ZeroMatrix(TDim, TDim - 1);
        noalias(jMaster) = ZeroMatrix(TDim, TDim - 1);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "MortarKinematicVariablesWithDerivatives with " << TNumNodes << " nodes (slave) and " << TNumNodesMaster << " nodes (master)";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MortarKinematicVariablesWithDerivatives with " << TNumNodes << " nodes (slave) and " << TNumNodesMaster << " nodes (master)" << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        // Print the base class data
        BaseClassType::PrintData(rOStream);

        // Print the additional variables
        rOStream << "DNDeMaster: " << DNDeMaster << "\n";
        rOStream << "DNDeSlave: " << DNDeSlave << "\n";
        rOStream << "jSlave: " << jSlave << "\n";
        rOStream << "jMaster: " << jMaster << std::endl;
    }

    ///@}
private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseClassType );
        rSerializer.save("DNDeMaster", DNDeMaster);
        rSerializer.save("DNDeSlave", DNDeSlave);
        rSerializer.save("jSlave", jSlave);
        rSerializer.save("jMaster", jMaster);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseClassType );
        rSerializer.load("DNDeMaster", DNDeMaster);
        rSerializer.load("DNDeSlave", DNDeSlave);
        rSerializer.load("jSlave", jSlave);
        rSerializer.load("jMaster", jMaster);
    }

    ///@}
}; // Class MortarKinematicVariablesWithDerivatives

/**
 * @class DerivativeData
 * @ingroup KratosCore
 * @brief This data will be used to compute the derivatives.
 * @details This class includes different information that is used in order to compute the directional derivatives in the mortar contact conditions
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template< const SizeType TDim, const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
class DerivativeData
{
public:
    ///@name Type Definitions
    ///@{

    // Geometry definition
    using GeometryType = Geometry<Node>;

    /// The definition of the index type
    using IndexType = std::size_t;

    // Auxiliary types
    using DummyBoundedMatrixType = BoundedMatrix<int, 1, 1>;

    using GeometryArraySlaveType = array_1d<double, TNumNodes>;

    using GeometryArrayMasterType = array_1d<double, TNumNodesMaster>;

    using GeometryDoFMatrixSlaveType = BoundedMatrix<double, TNumNodes, TDim>;

    using GeometryDoFMatrixMasterType = BoundedMatrix<double, TNumNodesMaster, TDim>;

    using GeometryMatrixType = BoundedMatrix<double, TNumNodes, TNumNodes>;

    using VertexDerivativesMatrixType = typename std::conditional<TNumNodes == 2, DummyBoundedMatrixType, BoundedMatrix<double, 3, 3>>::type;

    // Auxiliary sizes
    static const SizeType DummySize = 1;

    static const SizeType DoFSizeSlaveGeometry = (TNumNodes * TDim);

    static const SizeType DoFSizeMasterGeometry = (TNumNodesMaster * TDim);

    static const SizeType DoFSizePairedGeometry = DoFSizeSlaveGeometry + DoFSizeMasterGeometry;

    static const SizeType DoFSizeDerivativesDependence = (TDim == 2) ? DoFSizeSlaveGeometry : DoFSizePairedGeometry;

    static const SizeType DoFSizeDerivativesVertex = (TDim == 2) ? DummySize : DoFSizePairedGeometry;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    DerivativeData()= default;

    /// Destructor
    virtual ~DerivativeData()= default;

    ///@}
    ///@name Member Variables
    ///@{

    /// The penalty parameter
    array_1d<double, TNumNodes> PenaltyParameter;
    /// The scale factor
    double ScaleFactor;

    /// The normals of the nodes
    GeometryDoFMatrixSlaveType NormalSlave;
    GeometryDoFMatrixMasterType NormalMaster;

    /// Displacements and original coordinates
    GeometryDoFMatrixSlaveType X1, u1;
    GeometryDoFMatrixMasterType X2, u2;

    /// Jacobian derivatives
    array_1d<double, DoFSizeDerivativesDependence> DeltaDetjSlave;
    /// Dual shape function derivatives
    array_1d<GeometryArraySlaveType, DoFSizeDerivativesDependence> DeltaPhi;
    /// Shape function derivatives
    array_1d<GeometryArraySlaveType, DoFSizeDerivativesDependence> DeltaN1;
    array_1d<GeometryArrayMasterType, DoFSizeDerivativesDependence> DeltaN2;
    /// Normal derivatives
    array_1d<GeometryDoFMatrixSlaveType, DoFSizeSlaveGeometry> DeltaNormalSlave;
    array_1d<GeometryDoFMatrixMasterType, DoFSizeMasterGeometry> DeltaNormalMaster;
    /// Integration cell vertex derivatives
    array_1d<VertexDerivativesMatrixType, DoFSizeDerivativesVertex> DeltaCellVertex;

    /**
     * @brief Ae matrix.
     * @details This matrix corresponds with a coefficient matrix that relates the standard shape functions with the dual shape function
     * It is used to stabilize and for an assembling reduction
     * The name is taken from Popp's work
     */
    GeometryMatrixType Ae;

    /// Derivatives Ae
    array_1d<GeometryMatrixType, DoFSizeDerivativesDependence> DeltaAe;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Initializer method
     * @param SlaveGeometry The geometry of the slave
     * @param rCurrentProcessInfo The process info from the system
     */
    virtual void Initialize(
        const GeometryType& SlaveGeometry,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        // The normals of the nodes
        noalias(NormalSlave) = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry,  NORMAL, 0);

        // Displacements and velocities of the slave
        const IndexType step = (rCurrentProcessInfo[STEP] == 1) ? 0 : 1;
        noalias(u1) = step == 0 ?
                      MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry, DISPLACEMENT, 0) :
                      MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry, DISPLACEMENT, 0)
                    - MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry, DISPLACEMENT, 1);
        noalias(X1) = MortarUtilities::GetCoordinates<TDim,TNumNodes>(SlaveGeometry, false, step);

        // We get the ALM variables
        for (IndexType i = 0; i < TNumNodes; ++i)
            PenaltyParameter[i] = SlaveGeometry[i].GetValue(INITIAL_PENALTY);
        ScaleFactor = rCurrentProcessInfo[SCALE_FACTOR];

        // We initialize the derivatives
        const array_1d<double, TNumNodes> zero_vector_slave(TNumNodes, 0.0);
        const array_1d<double, TNumNodesMaster> zero_vector_master(TNumNodesMaster, 0.0);
        for (IndexType i = 0; i < DoFSizeSlaveGeometry; ++i) {
            DeltaDetjSlave[i] = 0.0;
            noalias(DeltaPhi[i]) = zero_vector_slave;
            noalias(DeltaN1[i]) = zero_vector_slave;
            noalias(DeltaN2[i]) = zero_vector_master;
        }
        if constexpr (TDim == 3) {
            for (IndexType i = 0; i < DoFSizeMasterGeometry; ++i) {
                DeltaDetjSlave[i + DoFSizeSlaveGeometry] = 0.0;
                noalias(DeltaPhi[i + DoFSizeSlaveGeometry]) = zero_vector_slave;
                noalias(DeltaN1[i + DoFSizeSlaveGeometry]) = zero_vector_slave;
                noalias(DeltaN2[i + DoFSizeSlaveGeometry]) = zero_vector_master;
            }
        }
    }

    /**
     * @brief This method reset to zero the cell vertex derivatives
     */
    virtual void ResetDerivatives()
    {
        // Derivatives
        if constexpr (TDim == 3) { // Derivative of the cell vertex
            // Auxiliary zero matrix
            const BoundedMatrix<double, 3, 3> aux_zero = ZeroMatrix(3, 3);
            for (IndexType i = 0; i < TNumNodes * TDim; ++i) {
                noalias(DeltaCellVertex[i]) = aux_zero;
            }
            for (IndexType i = 0; i < TNumNodesMaster * TDim; ++i) {
                noalias(DeltaCellVertex[i + TNumNodes * TDim]) = aux_zero;
            }
        }
    }

    /**
     * @brief Initialize the DeltaAe components
     */
    void InitializeDeltaAeComponents()
    {
        // Auxiliary zero matrix
        const GeometryMatrixType aux_zero = ZeroMatrix(TNumNodes, TNumNodes);

        // Ae
        noalias(Ae) = aux_zero;

        // Derivatives Ae
        for (IndexType i = 0; i < DoFSizeDerivativesDependence; ++i)
            noalias(DeltaAe[i]) = aux_zero;
    }

    /**
     * @brief Updating the Master pair
     * @param MasterGeometry The master geometry
     * @param rCurrentProcessInfo The process info from the system
     */
    virtual void UpdateMasterPair(
        const GeometryType& MasterGeometry,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        NormalMaster = MortarUtilities::GetVariableMatrix<TDim,TNumNodesMaster>(MasterGeometry,  NORMAL, 0);

        // Displacements, coordinates and normals of the master
        const IndexType step = (rCurrentProcessInfo[STEP] == 1) ? 0 : 1;
        noalias(u2) = step == 0 ?
                      MortarUtilities::GetVariableMatrix<TDim,TNumNodesMaster>(MasterGeometry, DISPLACEMENT, 0) :
                      MortarUtilities::GetVariableMatrix<TDim,TNumNodesMaster>(MasterGeometry, DISPLACEMENT, 0)
                    - MortarUtilities::GetVariableMatrix<TDim,TNumNodesMaster>(MasterGeometry, DISPLACEMENT, 1);
        noalias(X2) = MortarUtilities::GetCoordinates<TDim,TNumNodesMaster>(MasterGeometry, false, step);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "DerivativeData " << TDim << "D with " << TNumNodes << " nodes (slave) and " << TNumNodesMaster << " nodes (master)";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DerivativeData " << TDim << "D with " << TNumNodes << " nodes (slave) and " << TNumNodesMaster << " nodes (master)" << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        // Print the base class data
        rOStream << "PenaltyParameter: " << PenaltyParameter << "\n";
        rOStream << "ScaleFactor: " << ScaleFactor << "\n";
        rOStream << "NormalSlave: " << NormalSlave << "\n";
        rOStream << "NormalMaster: " << NormalMaster << "\n";
        rOStream << "X1: " << X1 << "\n";
        rOStream << "u1: " << u1 << "\n";
        rOStream << "X2: " << X2 << "\n";
        rOStream << "u2: " << u2 << "\n";
        for (IndexType i = 0; i < DoFSizeSlaveGeometry; ++i) {
            rOStream << "DeltaDetjSlave[" << i << "]: " << DeltaDetjSlave[i] << "\n";
            rOStream << "DeltaPhi[" << i << "]: " << DeltaPhi[i] << "\n";
            rOStream << "DeltaN1[" << i << "]: " << DeltaN1[i] << "\n";
            rOStream << "DeltaN2[" << i << "]: " << DeltaN2[i] << "\n";
            rOStream << "DeltaNormalSlave[" << i << "]: " << DeltaNormalSlave[i] << "\n";
            rOStream << "DeltaNormalMaster[" << i << "]: " << DeltaNormalMaster[i] << "\n";
            rOStream << "DeltaCellVertex[" << i << "]: " << DeltaCellVertex[i] << "\n";
        }
        rOStream << "Ae: "  << Ae  << "\n";
        for (IndexType i = 0; i < DoFSizeDerivativesDependence; ++i) {
            rOStream <<"DeltaAe["<< i<<"]: "<< DeltaAe[i]  << " ";
        }
        rOStream << std::endl;
    }

    ///@}
private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        rSerializer.save("PenaltyParameter", PenaltyParameter);
        rSerializer.save("ScaleFactor", ScaleFactor);
        rSerializer.save("NormalSlave", NormalSlave);
        rSerializer.save("NormalMaster", NormalMaster);
        rSerializer.save("X1", X1);
        rSerializer.save("u1", u1);
        rSerializer.save("X2", X2);
        rSerializer.save("u2", u2);
        rSerializer.save("DeltaDetjSlave", DeltaDetjSlave);
        rSerializer.save("DeltaPhi", DeltaPhi);
        rSerializer.save("DeltaN1", DeltaN1);
        rSerializer.save("DeltaN2", DeltaN2);
        rSerializer.save("DeltaNormalSlave", DeltaNormalSlave);
        rSerializer.save("DeltaNormalMaster", DeltaNormalMaster);
        rSerializer.save("DeltaCellVertex", DeltaCellVertex);
        rSerializer.save("Ae", Ae);
        rSerializer.save("DeltaAe", DeltaAe);
    }

    virtual void load(Serializer& rSerializer)
    {
        rSerializer.load("PenaltyParameter", PenaltyParameter);
        rSerializer.load("ScaleFactor", ScaleFactor);
        rSerializer.load("NormalSlave", NormalSlave);
        rSerializer.load("NormalMaster", NormalMaster);
        rSerializer.load("X1", X1);
        rSerializer.load("u1", u1);
        rSerializer.load("X2", X2);
        rSerializer.load("u2", u2);
        rSerializer.load("DeltaDetjSlave", DeltaDetjSlave);
        rSerializer.load("DeltaPhi", DeltaPhi);
        rSerializer.load("DeltaN1", DeltaN1);
        rSerializer.load("DeltaN2", DeltaN2);
        rSerializer.load("DeltaNormalSlave", DeltaNormalSlave);
        rSerializer.load("DeltaNormalMaster", DeltaNormalMaster);
        rSerializer.load("DeltaCellVertex", DeltaCellVertex);
        rSerializer.load("Ae", Ae);
        rSerializer.load("DeltaAe", DeltaAe);
    }

    ///@}
};  // Class DerivativeData

/**
 * @class DerivativeDataFrictional
 * @ingroup KratosCore
 * @brief This class is a derived class of DerivativeData
 * @details Includes additionally the derivatives necessary to compute the directional derivatives for the frictional conditions
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template< const SizeType TDim, const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
class DerivativeDataFrictional
    : public DerivativeData<TDim, TNumNodes, TNumNodesMaster>
{
public:
    ///@name Type Definitions
    ///@{

    // Geometry definition
    using GeometryType = Geometry<Node>;

    /// The base class type
    using BaseClassType = DerivativeData<TDim, TNumNodes, TNumNodesMaster>;

    /// The bounded matrix employed class
    using GeometryDoFMatrixSlaveType = BoundedMatrix<double, TNumNodes, TDim>;
    using GeometryDoFMatrixMasterType = BoundedMatrix<double, TNumNodesMaster, TDim>;

    // Size of DoFs of a not paired dependency
    static const SizeType DoFSizeSlaveGeometry = (TNumNodes * TDim);
    static const SizeType DoFSizeMasterGeometry = (TNumNodes * TDim);

    /// Size of DoFs of a paired dependency
    static const SizeType DoFSizePairedGeometry = DoFSizeSlaveGeometry + DoFSizeMasterGeometry;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    DerivativeDataFrictional()= default;

    /// Destructor
    ~DerivativeDataFrictional() override = default;

    ///@}
    ///@name Member Variables
    ///@{

    /// The ALM parameters
    double TangentFactor = 0.0;

    /// Displacements and velocities
    GeometryDoFMatrixSlaveType u1old;
    GeometryDoFMatrixMasterType u2old;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Initializer method
     * @param SlaveGeometry The geometry of the slave
     * @param rCurrentProcessInfo The process info from the system
     */
    void Initialize(
        const GeometryType& SlaveGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        BaseClassType::Initialize(SlaveGeometry, rCurrentProcessInfo);

        TangentFactor = rCurrentProcessInfo[TANGENT_FACTOR];

        noalias(u1old) = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry, DISPLACEMENT, 1) - MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry, DISPLACEMENT, 2);
    }

    /**
     * @brief Updating the Master pair
     * @param MasterGeometry The geometry of the master
     * @param rCurrentProcessInfo The process info from the system
     */
    void UpdateMasterPair(
        const GeometryType& MasterGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        BaseClassType::UpdateMasterPair(MasterGeometry, rCurrentProcessInfo);

        noalias(u2old) = MortarUtilities::GetVariableMatrix<TDim,TNumNodesMaster>(MasterGeometry, DISPLACEMENT, 1) - MortarUtilities::GetVariableMatrix<TDim,TNumNodesMaster>(MasterGeometry, DISPLACEMENT, 2);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "DerivativeDataFrictional " << TDim << "D with " << TNumNodes << " nodes (slave) and " << TNumNodesMaster << " nodes (master)";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "DerivativeDataFrictional " << TDim << "D with " << TNumNodes << " nodes (slave) and " << TNumNodesMaster << " nodes (master)" << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        // Print the base class data
        BaseClassType::PrintData(rOStream);

        // Print the additional variables
        rOStream << "TangentFactor: " << TangentFactor << "\n";
        rOStream << "u1old: " << u1old << "\n";
        rOStream << "u2old: " << u2old << std::endl;
    }

    ///@}
private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseClassType );
        rSerializer.save("TangentFactor", TangentFactor);
        rSerializer.save("u1old", u1old);
        rSerializer.save("u2old", u2old);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseClassType );
        rSerializer.load("TangentFactor", TangentFactor);
        rSerializer.load("u1old", u1old);
        rSerializer.load("u2old", u2old);
    }

    ///@}
};  // Class DerivativeDataFrictional

/**
 * @class MortarOperator
 * @ingroup KratosCore
 * @brief  This is the definition of the mortar operator according to the work of Alexander Popp: https://www.lnm.mw.tum.de/staff/alexander-popp/
 * @details In particular the thesis of contact mechanics based in mortar method available at: https://mediatum.ub.tum.de/?id=1109994
 * These mortar operator are assembled as mass matrices in order to transfer information from the slave side to the master side.
 * The operators are DOperator for the slave side and MOperator for master respectively.
 * In order to compute these operators, the shape functions of both domains are necessary (using the slave side as reference), as well as the integration weight and jacobian in the integration point.
 * Popp thesis pge 50 and following
 * @author Vicente Mataix Ferrandiz
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template< const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
class MortarOperator
{
public:
    ///@name Type Definitions
    ///@{

    /// The kinematic variables class
    using KinematicVariables = MortarKinematicVariables<TNumNodes, TNumNodesMaster>;

    /// The bounded matrix employed class
    using GeometryMatrixSlaveType = BoundedMatrix<double, TNumNodes, TNumNodes>;
    using GeometryMatrixMasterType = BoundedMatrix<double, TNumNodes, TNumNodesMaster>;

    /// Counted pointer of MortarOperator
    KRATOS_CLASS_POINTER_DEFINITION( MortarOperator );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    MortarOperator() = default;

    /// Destructor
    virtual ~MortarOperator() = default;

    ///@}
    ///@name Member Variables
    ///@{

    /// Mortar condition matrices - DOperator and MOperator
    GeometryMatrixSlaveType DOperator;
    GeometryMatrixMasterType MOperator;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method initialized the operators
     */
    void Initialize()
    {
        // We initialize the D and M operators
        noalias(DOperator) = ZeroMatrix(TNumNodes, TNumNodes);
        noalias(MOperator) = ZeroMatrix(TNumNodes, TNumNodesMaster);
    }

    /**
     * @brief It calculates the mortar operators. Popp thesis page 56, equation 3.31 and 3.32
     * @param rKinematicVariables Corresponds with the kinematic variables
     * @param rIntegrationWeight The corresponding integration weight
     */
    void CalculateMortarOperators(
        KinematicVariables& rKinematicVariables,
        const double rIntegrationWeight
        )
    {
        /* DEFINITIONS */
        const double det_j_slave = rKinematicVariables.DetjSlave;
        const Vector& phi_vector  = rKinematicVariables.PhiLagrangeMultipliers;
        const Vector& n1_vector   = rKinematicVariables.NSlave;
        const Vector& n2_vector   = rKinematicVariables.NMaster;

        for (IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave) {
            const double phi = phi_vector[i_slave];
            for (IndexType j_slave = 0; j_slave < TNumNodes; ++j_slave) {
                DOperator(i_slave, j_slave) += det_j_slave * rIntegrationWeight * phi * n1_vector[j_slave];
            }
            for (IndexType j_slave = 0; j_slave < TNumNodesMaster; ++j_slave) {
                MOperator(i_slave, j_slave) += det_j_slave * rIntegrationWeight * phi * n2_vector[j_slave];
            }
        }
    }

    /**
     * @brief It calculates the POperator (Inverse(D x M))
     * @details Popp thesis page 83 equation 3.88
     */
    GeometryMatrixMasterType ComputePOperator()
    {
        // We calculate the inverse of D operator
        double auxdet;
        GeometryMatrixSlaveType inv_D_operator;
        MathUtils<double>::InvertMatrix(DOperator, inv_D_operator, auxdet);

        // We calculate the P operator
        const GeometryMatrixMasterType POperator = prod(inv_D_operator, MOperator);

        return POperator;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "MortarOperator with " << TNumNodes << " nodes (slave) and " << TNumNodesMaster << " nodes (master)";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MortarOperator with " << TNumNodes << " nodes (slave) and " << TNumNodesMaster << " nodes (master)" << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        // Print the base class data
        rOStream << "DOperator: " << DOperator << "\n";
        rOStream << "MOperator: " << MOperator << std::endl;
    }

    ///@}
private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        rSerializer.save("DOperator", DOperator);
        rSerializer.save("MOperator", MOperator);
    }

    virtual void load(Serializer& rSerializer)
    {
        rSerializer.load("DOperator", DOperator);
        rSerializer.load("MOperator", MOperator);
    }

    ///@}
}; // Class MortarOperator

/**
 * @class MortarOperatorWithDerivatives
 * @ingroup KratosCore
 * @brief  This class derives from the MortarOperator class and it includes the derived operators.
 * @details The derived operators are defined in each DoF of each domain, which means TNumNodes x TDim x 2 derivatives definitions in order to compute all the necessary derivatives. Popp thesis page 102 and following
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TFrictional If the problem is frictional or not
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template< const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const SizeType TNumNodesMaster = TNumNodes>
class MortarOperatorWithDerivatives
    : public MortarOperator<TNumNodes, TNumNodesMaster>
{
public:
    ///@name Type Definitions
    ///@{

    /// The definition of the index type
    using IndexType = std::size_t;

    /// The base class type
    using BaseClassType = MortarOperator<TNumNodes, TNumNodesMaster>;

    /// The kinematic variables class
    using KinematicVariables = MortarKinematicVariables<TNumNodes, TNumNodesMaster>;

    /// The type for frictional derivative data
    using DerivativeDataFrictionalType = DerivativeDataFrictional<TDim, TNumNodes, TNumNodesMaster>;

    /// The type for frictionless derivative data
    using DerivativeFrictionalessDataType = DerivativeData<TDim, TNumNodes, TNumNodesMaster>;

    /// The type for derivative data, determined by the presence of friction
    using DerivativeDataType = typename std::conditional<TFrictional, DerivativeDataFrictionalType, DerivativeFrictionalessDataType>::type;

    /// The bounded matrix employed class for slave geometry
    using GeometryMatrixSlaveType = BoundedMatrix<double, TNumNodes, TNumNodes>;

    /// The bounded matrix employed class for master geometry
    using GeometryMatrixMasterType = BoundedMatrix<double, TNumNodes, TNumNodesMaster>;

    // Auxiliary sizes
    static const SizeType DoFSizeSlaveGeometry = (TNumNodes * TDim);
    static const SizeType DoFSizeMasterGeometry = (TNumNodesMaster * TDim);

    static const SizeType DoFSizePairedGeometry = DoFSizeSlaveGeometry + DoFSizeMasterGeometry;

    static const SizeType DoFSizeDerivativesDependence = (TDim == 2) ? DoFSizeSlaveGeometry : DoFSizePairedGeometry;

    /// Counted pointer of MortarOperatorWithDerivatives
    KRATOS_CLASS_POINTER_DEFINITION( MortarOperatorWithDerivatives );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    MortarOperatorWithDerivatives() = default;

    /// Destructor
    ~MortarOperatorWithDerivatives() override = default;

    ///@}
    ///@name Member Variables
    ///@{

    // D and M directional derivatives
    array_1d<GeometryMatrixSlaveType, DoFSizePairedGeometry> DeltaDOperator;
    array_1d<GeometryMatrixMasterType, DoFSizePairedGeometry> DeltaMOperator;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method initialized the operators
     */
    void Initialize()
    {
        BaseClassType::Initialize();

        // Auxiliary zero matrix
        const GeometryMatrixSlaveType aux_zero_slave = ZeroMatrix(TNumNodes, TNumNodes);
        const GeometryMatrixMasterType aux_zero_master = ZeroMatrix(TNumNodes, TNumNodesMaster);

        // We initialize the D and M derivatives operators
        for (IndexType i = 0; i < DoFSizeSlaveGeometry; ++i) {
            noalias(DeltaDOperator[i]) = aux_zero_slave;
            noalias(DeltaMOperator[i]) = aux_zero_master;
        }
        for (IndexType i = 0; i < DoFSizeMasterGeometry; ++i) {
            noalias(DeltaDOperator[i + DoFSizeSlaveGeometry]) = aux_zero_slave;
            noalias(DeltaMOperator[i + DoFSizeSlaveGeometry]) = aux_zero_master;
        }
    }

    /**
     * @brief It calculates the mortar operators
     * @details Popp thesis page 102 equation equation 4.32 and 4.33 / 4.37 and 4.38
     * @param rKinematicVariables Corresponds with the kinematic variables
     * @param rIntegrationWeight The corresponding integration weight
     */
    void CalculateDeltaMortarOperators(
        KinematicVariables& rKinematicVariables,
        DerivativeDataType& rDerivativeData,
        const double rIntegrationWeight
        )
    {
        /* DEFINITIONS */
        const double det_j_slave = rKinematicVariables.DetjSlave;
        const Vector& vector_phi = rKinematicVariables.PhiLagrangeMultipliers;
        const Vector& vector_n1  = rKinematicVariables.NSlave;
        const Vector& vector_n2  = rKinematicVariables.NMaster;

        // Derivatives
        const array_1d<double, DoFSizeDerivativesDependence>& delta_det_j_slave = rDerivativeData.DeltaDetjSlave;
        const array_1d<array_1d<double, TNumNodes >, DoFSizeDerivativesDependence>& delta_phi = rDerivativeData.DeltaPhi;
        const array_1d<array_1d<double, TNumNodes >, DoFSizeDerivativesDependence>& delta_n1  = rDerivativeData.DeltaN1;
        const array_1d<array_1d<double, TNumNodesMaster >, DoFSizeDerivativesDependence>& delta_n2  = rDerivativeData.DeltaN2;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const double phi = vector_phi[i_node];

            for (IndexType j_node = 0; j_node < TNumNodes; ++j_node) {
                const double n1 = vector_n1[j_node];

                BaseClassType::DOperator(i_node, j_node) += det_j_slave * rIntegrationWeight * phi * n1;

                for (IndexType i = 0; i < DoFSizeSlaveGeometry; ++i) {
                    DeltaDOperator[i](i_node, j_node) += delta_det_j_slave[i] * rIntegrationWeight * phi* n1
                                                       + det_j_slave * rIntegrationWeight * delta_phi[i][i_node] * n1
                                                       + det_j_slave * rIntegrationWeight * phi* delta_n1[i][j_node];
                }
                if constexpr (TDim == 3) {
                    for (IndexType i = DoFSizeSlaveGeometry; i < DoFSizePairedGeometry; ++i) {
                        DeltaDOperator[i](i_node, j_node) += det_j_slave * rIntegrationWeight * phi * delta_n1[i][j_node];
                        DeltaDOperator[i](i_node, j_node) += delta_det_j_slave[i] * rIntegrationWeight * phi * n1;
                        DeltaDOperator[i](i_node, j_node) += det_j_slave * rIntegrationWeight * delta_phi[i][i_node] * n1;
                    }
                }
            }
            for (IndexType j_node = 0; j_node < TNumNodesMaster; ++j_node) {
                const double n2 = vector_n2[j_node];

                BaseClassType::MOperator(i_node, j_node) += det_j_slave * rIntegrationWeight * phi * n2;

                for (IndexType i = 0; i < DoFSizeSlaveGeometry; ++i) {
                    DeltaMOperator[i](i_node, j_node) += delta_det_j_slave[i] * rIntegrationWeight * phi* n2
                                                       + det_j_slave * rIntegrationWeight * delta_phi[i][i_node] * n2
                                                       + det_j_slave * rIntegrationWeight * phi* delta_n2[i][j_node];
                }
                if constexpr (TDim == 3) {
                    for (IndexType i = DoFSizeSlaveGeometry; i < DoFSizePairedGeometry; ++i) {
                        DeltaMOperator[i](i_node, j_node) += det_j_slave * rIntegrationWeight * phi * delta_n2[i][j_node];
                        DeltaMOperator[i](i_node, j_node) += delta_det_j_slave[i] * rIntegrationWeight * phi * n2;
                        DeltaMOperator[i](i_node, j_node) += det_j_slave * rIntegrationWeight * delta_phi[i][i_node] * n2;
                    }
                }
            }
        }
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "MortarOperatorWithDerivatives in " << TDim << "D with " << TNumNodes << " nodes (slave) and " << TNumNodesMaster << " nodes (master)";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MortarOperatorWithDerivatives in " << TDim << "D with " << TNumNodes << " nodes (slave) and " << TNumNodesMaster << " nodes (master)" << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        // Print the base class data
        BaseClassType::PrintData(rOStream);

        // Print the additional variables
        for (IndexType i = 0; i < TNumNodes * TDim; ++i) {
            rOStream << "DeltaDOperator_" << i << ": " << DeltaDOperator[i] << std::endl;
            rOStream << "DeltaMOperator_" << i << ": " << DeltaMOperator[i] << std::endl;
        }
    }

    ///@}
private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseClassType );
        rSerializer.save("DeltaDOperator", DeltaDOperator);
        rSerializer.save("DeltaMOperator", DeltaMOperator);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseClassType );
        rSerializer.load("DeltaDOperator", DeltaDOperator);
        rSerializer.load("DeltaMOperator", DeltaMOperator);
    }

    ///@}
}; // Class MortarOperatorWithDerivatives

/**
 * @class DualLagrangeMultiplierOperators
 * @ingroup KratosCore
 * @brief  This is the definition dual lagrange multiplier operators according to the work of Alexander Popp: https://www.lnm.mw.tum.de/staff/alexander-popp/
 * @details In particular the thesis of contact mechanics based in mortar method available at: https://mediatum.ub.tum.de/?id=1109994
 * In order to compute the dual LM shape function the Ae operator must be computed, which depends of the Me and De operators. Phi = Ae * NSlave.  In a similar way to the mortar operators, the De corresponds with a diagonal operator and Me with a sparse operator respectively. Ae = De * inv(Me)
 * Popp thesis page 69 and following
 * @author Vicente Mataix Ferrandiz
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template< const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
class DualLagrangeMultiplierOperators
{
public:
    ///@name Type Definitions
    ///@{

    /// The definition of the index type
    using IndexType = std::size_t;

    /// The kinematic variables class
    using KinematicVariables = MortarKinematicVariables<TNumNodes, TNumNodesMaster>;

    /// The bounded matrix employed class
    using GeometryMatrixType = BoundedMatrix<double, TNumNodes, TNumNodes>;

    /// Counted pointer of DualLagrangeMultiplierOperators
    KRATOS_CLASS_POINTER_DEFINITION( DualLagrangeMultiplierOperators );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    DualLagrangeMultiplierOperators() = default;

    /// Destructor
    virtual ~DualLagrangeMultiplierOperators() = default;

    ///@}
    ///@name Member Variables
    ///@{

    /// The auxiliary operators needed to build the dual LM Ae operator
    GeometryMatrixType Me, De;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method initialized the operators
     */
    void Initialize()
    {
        // We initialize the De and Me operators
        noalias(Me) = ZeroMatrix(TNumNodes, TNumNodes);
        noalias(De) = ZeroMatrix(TNumNodes, TNumNodes);
    }

    /**
     * @brief Calculates the Ae components necessary to compute the Phi_LagrangeMultipliers shape functions.
     * @details For that it integrates De and Me. Popp thesis page 70 eq. 3.65
     * @param rKinematicVariables The kinematic variables
     * @param rIntegrationWeight The integration weight considered
     */
    void CalculateAeComponents(
        KinematicVariables& rKinematicVariables,
        const double rIntegrationWeight
        )
    {
        /* DEFINITIONS */
        const Vector& n1 = rKinematicVariables.NSlave;
        const double det_j = rKinematicVariables.DetjSlave;

        noalias(De) += rIntegrationWeight * (ComputeDe(n1, det_j));
        noalias(Me) += rIntegrationWeight * det_j * outer_prod(n1, n1);
    }

    /**
     * @brief Calculates the matrix Ae. To avoid problems in the inversion the matrix is normalized
     * @details Popp thesis page 70. Equation 3.65
     * @param Ae The dual Lagrange Multiplier operator
     */
    bool CalculateAe(GeometryMatrixType& Ae)
    {
        const double tolerance = std::numeric_limits<double>::epsilon();

        // We compute the norm
        const double norm_me = norm_frobenius(Me);

        // Now we normalize the matrix
        if (norm_me >= tolerance) {
            const GeometryMatrixType normalized_Me = Me/norm_me;

            // We compute the normalized inverse
            double aux_det;
            GeometryMatrixType normalized_inv_Me;
            MathUtils<double>::InvertMatrix(normalized_Me, normalized_inv_Me, aux_det, -1.0);
            const bool good_condition_number = MathUtils<double>::CheckConditionNumber(normalized_Me, normalized_inv_Me, tolerance, false);
            if (good_condition_number) {
                noalias(Ae) = (1.0/norm_me) * prod(De, normalized_inv_Me);
                return true;
            }
        }

        noalias(Ae) = IdentityMatrix(TNumNodes);
        return false;
    }

    /**
     * @brief Calculates the matrix De
     * @param N1 The shape function
     * @param detJ The jacobian of the geometry
     */
    template<class TArray>
    GeometryMatrixType ComputeDe(
        const TArray& N1,
        const double detJ
        ) const
    {
        GeometryMatrixType De;

        for (IndexType i = 0; i < TNumNodes; ++i) {
            for (IndexType j = 0; j < TNumNodes; ++j) {
                if (i == j) De(i,i) = detJ * N1[i];
                else De(i,j) = 0.0;
            }
        }

        return De;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "DualLagrangeMultiplierOperators with " << TNumNodes << " nodes (slave) and " << TNumNodesMaster << " nodes (master)";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DualLagrangeMultiplierOperators with " << TNumNodes << " nodes (slave) and " << TNumNodesMaster << " nodes (master)" << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const 
    {
        // Print the base class data
        rOStream << "Me: " << Me << std::endl;
        rOStream << "De: " << De << std::endl;
    }

    ///@}
private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        rSerializer.save("Me", Me);
        rSerializer.save("De", De);
    }

    virtual void load(Serializer& rSerializer)
    {
        rSerializer.load("Me", Me);
        rSerializer.load("De", De);
    }

    ///@}

}; // Class DualLagrangeMultiplierOperators

/**
 * @class DualLagrangeMultiplierOperatorsWithDerivatives
 * @ingroup KratosCore
 * @brief  This is the definition dual lagrange multiplier operators including the derivatives.
 * @details It is based in the same work as the previous class. In this case it computes the derivatives in order to compute the directionald erivative of the dual shape functions. Popp thesis page 111 and following
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TFrictional If the problem is frictional or not
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template< const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const SizeType TNumNodesMaster = TNumNodes>
class DualLagrangeMultiplierOperatorsWithDerivatives
    : public DualLagrangeMultiplierOperators<TNumNodes, TNumNodesMaster>
{
public:
    ///@name Type Definitions
    ///@{

    /// The definition of the index type
    using IndexType = std::size_t;

    /// The base class type
    using BaseClassType = DualLagrangeMultiplierOperators<TNumNodes, TNumNodesMaster>;

    /// The kinematic variables type with derivatives
    using KinematicVariablesType = MortarKinematicVariablesWithDerivatives<TDim, TNumNodes, TNumNodesMaster>;

    /// The type for frictional derivative data
    using DerivativeDataFrictionalType = DerivativeDataFrictional<TDim, TNumNodes, TNumNodesMaster>;

    /// The type for frictionless derivative data
    using DerivativeFrictionalessDataType = DerivativeData<TDim, TNumNodes, TNumNodesMaster>;

    /// The type for derivative data, determined by the presence of friction
    using DerivativeDataType = typename std::conditional<TFrictional, DerivativeDataFrictionalType, DerivativeFrictionalessDataType>::type;

    /// The bounded matrix employed class for geometry
    using GeometryMatrixType = BoundedMatrix<double, TNumNodes, TNumNodes>;

    // Auxiliary sizes
    static const SizeType DoFSizeSlaveGeometry = (TNumNodes * TDim);

    static const SizeType DoFSizeMasterGeometry = (TNumNodesMaster * TDim);

    static const SizeType DoFSizePairedGeometry = DoFSizeSlaveGeometry + DoFSizeMasterGeometry;

    static const SizeType DoFSizeDerivativesDependence = (TDim == 2) ? DoFSizeSlaveGeometry : DoFSizePairedGeometry;

    /// Counted pointer of DualLagrangeMultiplierOperatorsWithDerivatives
    KRATOS_CLASS_POINTER_DEFINITION( DualLagrangeMultiplierOperatorsWithDerivatives );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    DualLagrangeMultiplierOperatorsWithDerivatives() = default;

    /// Destructor
    ~DualLagrangeMultiplierOperatorsWithDerivatives() override = default;

    ///@}
    ///@name Member Variables
    ///@{

    // Derivatives matrices
    array_1d<GeometryMatrixType, DoFSizeDerivativesDependence> DeltaMe, DeltaDe;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method initialized the operators
     */
    void Initialize()
    {
        BaseClassType::Initialize();

        // Derivatives matrices
	    const BoundedMatrix<double, TNumNodes, TNumNodes> zeromatrix = ZeroMatrix(TNumNodes, TNumNodes);
        for (IndexType i = 0; i < DoFSizeDerivativesDependence; ++i) {
            noalias(DeltaMe[i]) = zeromatrix;
            noalias(DeltaDe[i]) = zeromatrix;
        }
    }

    /**
     * @brief Calculates the Ae components and its derivatives necessary to compute the Phi_LagrangeMultipliers shape functions.
     * @details Popp thesis page 112 eq. 4.59
     * @param rKinematicVariables The kinematic variables
     * @param rDerivativeData The data containing the derivatives
     * @param rIntegrationWeight The integration weight considered
     */
    void CalculateDeltaAeComponents(
        KinematicVariablesType& rKinematicVariables,
        DerivativeDataType& rDerivativeData,
        const double rIntegrationWeight
        )
    {
        /* DEFINITIONS */
        const double det_j_slave = rKinematicVariables.DetjSlave;
        const Vector& n1 = rKinematicVariables.NSlave;

        BaseClassType::CalculateAeComponents(rKinematicVariables, rIntegrationWeight);

        for (IndexType i = 0; i < DoFSizeDerivativesDependence; ++i) {
            const double delta_det_j = rDerivativeData.DeltaDetjSlave[i];
            const array_1d<double, TNumNodes>& delta_n1 = rDerivativeData.DeltaN1[i];

            noalias(DeltaDe[i]) += rIntegrationWeight * this->ComputeDe( n1, delta_det_j )
                                +  rIntegrationWeight * this->ComputeDe( delta_n1, det_j_slave );

            noalias(DeltaMe[i]) += rIntegrationWeight * delta_det_j * outer_prod(n1, n1)
                                +  rIntegrationWeight * det_j_slave * (outer_prod(delta_n1, n1) + outer_prod(n1, delta_n1));
        }
    }

    /**
     * @brief Calculates the matrix DeltaAe.
     * @details Popp thesis page 112 equation 4.58
     * @param rDerivativeData The data containing the derivatives
     */
    bool CalculateDeltaAe(DerivativeDataType& rDerivativeData)
    {
        const double tolerance = std::numeric_limits<double>::epsilon();

        // We compute the norm
        const double norm_Me = norm_frobenius(BaseClassType::Me);

        // Now we normalize the matrix
        const GeometryMatrixType normalized_Me = BaseClassType::Me/norm_Me;

        // We compute the normalized inverse
        double aux_det;
        GeometryMatrixType normalized_inv_Me;
        MathUtils<double>::InvertMatrix(normalized_Me, normalized_inv_Me, aux_det, -1.0);
        const bool good_condition_number = MathUtils<double>::CheckConditionNumber(normalized_Me, normalized_inv_Me, tolerance, false);
        if (!good_condition_number) {
            return false;
        }

        // Now we compute the inverse
        const GeometryMatrixType inv_Me = normalized_inv_Me/norm_Me;

        noalias(rDerivativeData.Ae) = prod(BaseClassType::De, inv_Me);

        array_1d<GeometryMatrixType , DoFSizeDerivativesDependence>& delta_Ae = rDerivativeData.DeltaAe;

        for (IndexType i = 0; i < DoFSizeDerivativesDependence; ++i) {
            const GeometryMatrixType aux_matrix = DeltaDe[i] - prod(rDerivativeData.Ae, DeltaMe[i]);
            noalias(delta_Ae[i]) = prod(aux_matrix, inv_Me);
        }

        return true;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "DualLagrangeMultiplierOperatorsWithDerivatives in " << TDim << "D with " << TNumNodes << " nodes (slave) and " << TNumNodesMaster << " nodes (master). Frictional: " << TFrictional;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "DualLagrangeMultiplierOperatorsWithDerivatives in " << TDim << "D with " << TNumNodes << " nodes (slave) and " << TNumNodesMaster << " nodes (master). Frictional: " << TFrictional << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        // Print the base class data
        BaseClassType::PrintData(rOStream);

        // Derivatives matrices
        for (IndexType i = 0; i < DoFSizeDerivativesDependence; ++i) {
            rOStream << "DeltaMe_" << i << ": " << DeltaMe[i] << std::endl;
            rOStream << "DeltaDe_" << i << ": " << DeltaDe[i] << std::endl;
        }
    }

    ///@}
private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseClassType );
        rSerializer.save("DeltaMe", DeltaMe);
        rSerializer.save("DeltaDe", DeltaDe);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseClassType );
        rSerializer.load("DeltaMe", DeltaMe);
        rSerializer.load("DeltaDe", DeltaDe);
    }

    ///@}

}; // Class DualLagrangeMultiplierOperatorsWithDerivatives

/**
 * @class PointBelong
 * @ingroup KratosCore
 * @brief Custom Point container to be used by the mapper
 * @details This point which is a derived class of the standard point, contains the variable mBelongs. This variable is a "hash" that can be used to determine where in which intersections the point belongs
 * @author Vicente Mataix Ferrandiz
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template<const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
class PointBelong
    : public Point
{
public:
    ///@name Type Definitions
    ///@{

    /// The belonging type
    using BelongType = typename std::conditional<TNumNodes == 2, PointBelongsLine2D2N, typename std::conditional<TNumNodes == 3, typename std::conditional<TNumNodesMaster == 3, PointBelongsTriangle3D3N, PointBelongsTriangle3D3NQuadrilateral3D4N>::type, typename std::conditional<TNumNodesMaster == 3, PointBelongsQuadrilateral3D4NTriangle3D3N, PointBelongsQuadrilateral3D4N>::type>::type>::type;

    /// Counted pointer of PointBelong
    KRATOS_CLASS_POINTER_DEFINITION( PointBelong );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    PointBelong():
        Point()
    {}

    /// Constructor with coordinates
    PointBelong(const array_1d<double, 3> Coords):
        Point(Coords)
    {}

    /// Constructor with coordinates and belonging type
    PointBelong(const array_1d<double, 3> Coords, const BelongType& ThisBelongs):
        Point(Coords),
        mBelongs(ThisBelongs)
    {}

    /// Destructor.
    ~PointBelong() override= default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method allows to set where the point belongs
     */
    void SetBelong(BelongType ThisBelongs)
    {
        mBelongs = ThisBelongs;
    }

    /**
     * @brief This method recovers where the point belongs
     */
    BelongType GetBelong() const
    {
        return mBelongs;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "PointBelong with " << TNumNodes << " nodes (slave) and " << TNumNodesMaster << " nodes (master)";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PointBelong with " << TNumNodes << " nodes (slave) and " << TNumNodesMaster << " nodes (master)" << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        // Print the base class data
        Point::PrintData(rOStream);

        // Print the belonging type
        rOStream << "BelongType:" << static_cast<int>(mBelongs) << std::endl;
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    BelongType mBelongs; /// To know if the point belongs to the master/slave/intersection (just 3D) side

    ///@}
}; // Class PointBelong

///@}
///@name Input and output
///@{

// MortarKinematicVariables class

/// input stream function
template< const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
inline std::istream& operator >> (std::istream& rIStream,
                                  MortarKinematicVariables<TNumNodes,TNumNodesMaster>& rThis)
                                  {return rIStream;};

/// output stream function
template< const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MortarKinematicVariables<TNumNodes,TNumNodesMaster>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

// MortarKinematicVariablesWithDerivatives class

/// input stream function
template< const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
inline std::istream& operator >> (std::istream& rIStream,
                                  MortarKinematicVariablesWithDerivatives<TNumNodes,TNumNodesMaster>& rThis)
                                  {return rIStream;};

/// output stream function
template< const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MortarKinematicVariablesWithDerivatives<TNumNodes,TNumNodesMaster>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

// DerivativeData class

/// input stream function
template< const SizeType TDim, const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
inline std::istream& operator >> (std::istream& rIStream,
                                  DerivativeData<TDim,TNumNodes,TNumNodesMaster>& rThis)
                                  {return rIStream;};

/// output stream function
template< const SizeType TDim, const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DerivativeData<TDim,TNumNodes,TNumNodesMaster>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

// DerivativeDataFrictional class

/// input stream function
template< const SizeType TDim, const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
inline std::istream& operator >> (std::istream& rIStream,
                                  DerivativeDataFrictional<TDim,TNumNodes,TNumNodesMaster>& rThis)
                                  {return rIStream;};

/// output stream function
template< const SizeType TDim, const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DerivativeDataFrictional<TDim,TNumNodes,TNumNodesMaster>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

// MortarOperator class

/// input stream function
template< const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
inline std::istream& operator >> (std::istream& rIStream,
                                  MortarOperator<TNumNodes,TNumNodesMaster>& rThis)
                                  {return rIStream;};

/// output stream function
template< const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MortarOperator<TNumNodes,TNumNodesMaster>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

// MortarOperatorWithDerivatives class

/// input stream function
template< const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const SizeType TNumNodesMaster = TNumNodes>
inline std::istream& operator >> (std::istream& rIStream,
                                  MortarOperatorWithDerivatives<TDim,TNumNodes,TFrictional,TNumNodesMaster>& rThis)
                                  {return rIStream;};

/// output stream function
template< const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const SizeType TNumNodesMaster = TNumNodes>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MortarOperatorWithDerivatives<TDim,TNumNodes,TFrictional,TNumNodesMaster>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

// DualLagrangeMultiplierOperators class

/// input stream function
template< const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
inline std::istream& operator >> (std::istream& rIStream,
                                  DualLagrangeMultiplierOperators<TNumNodes,TNumNodesMaster>& rThis)
                                  {return rIStream;};

/// output stream function
template< const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DualLagrangeMultiplierOperators<TNumNodes,TNumNodesMaster>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

// DualLagrangeMultiplierOperatorsWithDerivatives class

/// input stream function
template< const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const SizeType TNumNodesMaster = TNumNodes>
inline std::istream& operator >> (std::istream& rIStream,
                                  DualLagrangeMultiplierOperatorsWithDerivatives<TDim,TNumNodes,TFrictional,TNumNodesMaster>& rThis)
                                  {return rIStream;};

/// output stream function
template< const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const SizeType TNumNodesMaster = TNumNodes>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DualLagrangeMultiplierOperatorsWithDerivatives  <TDim,TNumNodes,TFrictional,TNumNodesMaster>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

// PointBelong class

/// input stream function
template< const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
inline std::istream& operator >> (std::istream& rIStream,
                                  PointBelong<TNumNodes,TNumNodesMaster>& rThis)
                                  {return rIStream;};

/// output stream function
template< const SizeType TNumNodes, const SizeType TNumNodesMaster = TNumNodes>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const PointBelong<TNumNodes,TNumNodesMaster>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}// namespace Kratos.