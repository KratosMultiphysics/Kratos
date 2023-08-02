//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include <filesystem>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/io.h"

namespace Kratos
{
///@addtogroup ApplicationNameApplication
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

/// This class reads from STL file format and creates triangular elements in given model_part
/** The current version only reads triangles from the STL and not higher order polygons
 * The nodes corresponging to given vertices are not collapsed
 * A SubModelPart for each additional solid block will be created
 * For definition STL format please check https://en.wikipedia.org/wiki/STL_(file_format)
 * A sample file format with 3 triangles:
   solid 3 triangles
   facet normal  1.000000 0.000000 0.000000 
       outer loop 
          vertex 0.1 -2.56114e-08 0.1
          vertex 0.1 -0.499156 -0.0352136
          vertex 0.1 -0.473406 -0.0446259
       endloop 
   endfacet 
   facet normal  1.000000 -0.000000 0.000000 
       outer loop 
          vertex 0.1 -0.473406 -0.0446259
          vertex 0.1 -0.447464 -0.0534931
          vertex 0.1 -2.56114e-08 0.1
       endloop 
   endfacet 
   facet normal  1.000000 0.000000 0.000000 
       outer loop 
          vertex 0.1 -0.6 0.1
          vertex 0.1 -0.524702 -0.0252604
          vertex 0.1 -0.499156 -0.0352136
       endloop 
   endfacet 
   endsolid 3 triangles

*/
class KRATOS_API(KRATOS_CORE) StlIO 
    : public IO
{
public:
    ///@name Type Definitions
    ///@{

    using GeometriesMapType = ModelPart::GeometriesMapType;
    using NodesArrayType = Element::NodesArrayType;

    /// Pointer definition of StlIO
    KRATOS_CLASS_POINTER_DEFINITION(StlIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with filename, and open option with default being read mode
    StlIO(
        std::filesystem::path const& Filename,
        Parameters ThisParameters = Parameters());

    /// Constructor with stream.
    StlIO(
        Kratos::shared_ptr<std::iostream> pInputStream,
        Parameters ThisParameters = Parameters());

    /// Destructor.
    virtual ~StlIO(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    static Parameters GetDefaultParameters();

    void ReadModelPart(ModelPart & rThisModelPart) override;

    void WriteModelPart(const ModelPart & rThisModelPart) override;

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

    Parameters mParameters;

    std::size_t mNextNodeId = 0;
    std::size_t mNextElementId = 0;
    std::size_t mNextConditionId = 0;

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

    Kratos::shared_ptr<std::iostream> mpInputStream;
    Flags mOptions;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void ReadSolid(
        ModelPart & rThisModelPart,
        const std::function<void(ModelPart&, NodesArrayType&)>& rCreateEntityFunctor );

    void ReadFacet(
        ModelPart & rThisModelPart,
        const std::function<void(ModelPart&, NodesArrayType&)>& rCreateEntityFunctor);

    void ReadLoop(
        ModelPart & rThisModelPart,
        const std::function<void(ModelPart&, NodesArrayType&)>& rCreateEntityFunctor);

    Point ReadPoint();

    void ReadKeyword(std::string const& Keyword);

    template<class TContainerType>
    void WriteEntityBlock(const TContainerType& rThisEntities);

    void WriteGeometryBlock(const GeometriesMapType& rThisGeometries);

    void WriteFacet(const GeometryType & rGeom);

    bool IsValidGeometry(
        const Geometry<Node>& rGeometry,
        std::size_t& rNumDegenerateGeos) const;

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
    StlIO& operator=(StlIO const& rOther);

    /// Copy constructor.
    StlIO(StlIO const& rOther);

    ///@}

}; // Class StlIO

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                StlIO& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const StlIO& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.
