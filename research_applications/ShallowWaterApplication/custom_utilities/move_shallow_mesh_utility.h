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

#ifndef KRATOS_MOVE_MESH_UTILITY_H_INCLUDED
#define KRATOS_MOVE_MESH_UTILITY_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/binbased_fast_point_locator.h"


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

/**
 * @class MoveShallowMeshUtility
 * @ingroup KratosShallowWaterApplication
 * @brief Tools for lagrangian computations
 * @details Move the computational mesh over a background mesh and map data between them
 * @author Miguel Maso Sotomayor
 */
class KRATOS_API(SHALLOW_WATER_APPLICATION) MoveShallowMeshUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MoveShallowMeshUtility
    KRATOS_CLASS_POINTER_DEFINITION(MoveShallowMeshUtility);

    typedef Node NodeType;

    typedef Geometry<NodeType> GeometryType;

    typedef BinBasedFastPointLocator<2>::ResultIteratorType ResultIteratorType;

    typedef BinBasedFastPointLocator<2>::ResultContainerType ResultContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * @brief Constructor
    */
    MoveShallowMeshUtility(
        ModelPart& rLagrangianModelPart,
        ModelPart& rEulerianModelPart,
        Parameters ThisParameters);

    /**
    * @brief Destructor
    */
    virtual ~MoveShallowMeshUtility(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    int Check();

    void Initialize();

    void MoveMesh();

    void MapResults();

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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "MoveShallowMeshUtility";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << Info();}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


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

    ModelPart& mrLagrangianModelPart;
    ModelPart& mrEulerianModelPart;

    BinBasedFastPointLocator<2> mLagrangianSearchStructure;
    BinBasedFastPointLocator<2> mEulerianSearchStructure;
    int mMaxResults;

    std::vector<const Variable<double>*> mScalarVariablesToLagrangian;
    std::vector<const Variable<array_1d<double,3>>*> mVectorVariablesToLagrangian;

    std::vector<const Variable<double>*> mScalarVariablesToEulerian;
    std::vector<const Variable<array_1d<double,3>>*> mVectorVariablesToEulerian;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    const Parameters GetDefaultParameters() const;

    template<class TDataType>
    void FillVariablesList(std::vector<const Variable<TDataType>*>& rVariablesList, const Parameters VariablesNames);

    bool MoveNode(
        NodeType& rNode,
        double Dt,
        Vector& rN,
        Element::Pointer& pElement,
        ResultIteratorType& rResultBegin);

    void MapToLagrangian(
        NodeType& rNode,
        const Vector& rN,
        const Element::Pointer pElement);

    void MapToEulerian(
        NodeType& rNode,
        const Vector& rN,
        const Element::Pointer pElement,
        const bool IsFound);

    template<class TDataType>
    void InterpolateVariable(
        NodeType& rNode,
        const Vector& rN,
        const GeometryType& rGeometry,
        const Variable<TDataType>& rVariable);

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    MoveShallowMeshUtility& operator=(MoveShallowMeshUtility const& rOther);

    /// Copy constructor.
    MoveShallowMeshUtility(MoveShallowMeshUtility const& rOther);

    ///@}

}; // Class MoveShallowMeshUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                MoveShallowMeshUtility& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const MoveShallowMeshUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MOVE_MESH_UTILITY_H_INCLUDED  defined
