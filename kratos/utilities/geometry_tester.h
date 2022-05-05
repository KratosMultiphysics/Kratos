//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_GEOMETRY_TEST_H_INCLUDED )
#define  KRATOS_GEOMETRY_TEST_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/model.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@}
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
 * @class GeometryTesterUtility
 * @ingroup KratosCore
 * @brief This utility tests the geometries
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
*/
class GeometryTesterUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometryTesterUtility
    KRATOS_CLASS_POINTER_DEFINITION(GeometryTesterUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    GeometryTesterUtility() {}

    /// Destructor.
    virtual ~GeometryTesterUtility() {};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    bool RunTest(Model& rModel);

    bool TestTetrahedra3D4N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    bool TestTetrahedra3D10N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    bool TestTriangle2D3N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    bool TestTriangle2D6N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    bool TestQuadrilateral2D4N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    bool TestQuadrilateral2D9N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    bool TestQuadrilateralInterface2D4N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    bool TestHexahedra3D8N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    bool TestHexahedra3D20N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    bool TestHexahedra3D27N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    bool TestHexahedraInterface3D8N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    bool TestPrism3D6N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

   bool TestPrism3D15N(
       ModelPart& rModelPart,
       std::stringstream& rErrorMessage
       );

    bool TestPrismInterface3D6N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

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
        buffer << "GeometryTesterUtility" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GeometryTesterUtility";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

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

    void GenerateNodes(ModelPart& rModelPart);

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
    bool VerifyAreaByIntegration( 
        Geometry<Node<3> >& geom, 
        Geometry<Node<3> >::IntegrationMethod ThisMethod, 
        const double reference_area, 
        std::stringstream& rErrorMessage
        );

    //here we verify that a  "displacement field" which varies linearly in space, produces the expected strain distribution.
    //this shall be considered a test for shape function derivatives
    void VerifyStrainExactness( 
        Geometry<Node<3> >& geom,  
        Geometry<Node<3> >::IntegrationMethod ThisMethod, 
        std::stringstream& rErrorMessage
        );

    //this function computes the linear strain matrix - useful to verify that a constant strain can be correctly reproduced
    void CalculateB(
        Matrix& B,
        Matrix& DN_DX,
        const std::size_t number_of_nodes,
        const std::size_t dimension
        );

    bool VerifyIsInside(
        Geometry< Node<3> >& geom,
        Geometry< Node<3> >::CoordinatesArrayType& global_coordinates,
        bool expected_result,
        std::stringstream& rErrorMessage
        );

    bool VerfiyShapeFunctionsValues(
        Geometry< Node<3> >& geom,
        Geometry< Node<3> >::CoordinatesArrayType& global_coordinates,
        std::stringstream& rErrorMessage
        );

    std::string GetIntegrationName(Geometry< Node<3> >& geom, Geometry<Node<3> >::IntegrationMethod ThisMethod);

    std::string GetGeometryName(Geometry< Node<3> >& geom);

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
    GeometryTesterUtility& operator=(GeometryTesterUtility const& rOther) {return *this;}

    /// Copy constructor.
    GeometryTesterUtility(GeometryTesterUtility const& rOther) {}


    ///@}

}; // Class GeometryTesterUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  GeometryTesterUtility& rThis) {
  return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GeometryTesterUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_TEST_H_INCLUDED  defined
