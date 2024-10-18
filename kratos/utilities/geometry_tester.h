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
class KRATOS_API(KRATOS_CORE) GeometryTesterUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometryTesterUtility
    KRATOS_CLASS_POINTER_DEFINITION(GeometryTesterUtility);

    /// Node type
    typedef Node NodeType;

    /// Geometry type
    typedef Geometry<NodeType> GeometryType;

    typedef typename GeometryType::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;

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

    /**
     * @brief This function tests the geometries
     * @param rModel A model containing a model part
     * @return true If teh test fails, true otherwise¡
     */
    bool RunTest(Model& rModel);

    /**
     * @brief This function tests the Tetrahedra3D4N
     * @param rModelPart Model part containing nodes
     * @return true If teh test fails, true otherwise¡
     */
    bool TestTetrahedra3D4N(ModelPart& rModelPart) {std::stringstream ss; return StreamTestTetrahedra3D4N(rModelPart, ss);}

    /**
     * @brief This function tests the Tetrahedra3D4N
     * @param rModelPart Model part containing nodes
     * @param rErrorMessage The error message
     * @return true If teh test fails, true otherwise¡
     */
    bool StreamTestTetrahedra3D4N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    /**
     * @brief This function tests the Tetrahedra3D10N
     * @param rModelPart Model part containing nodes
     * @return true If teh test fails, true otherwise¡
     */
    bool TestTetrahedra3D10N(ModelPart& rModelPart) {std::stringstream ss; return StreamTestTetrahedra3D10N(rModelPart, ss);}

    /**
     * @brief This function tests the Tetrahedra3D10N
     * @param rModelPart Model part containing nodes
     * @param rErrorMessage The error message
     * @return true If teh test fails, true otherwise¡
     */
    bool StreamTestTetrahedra3D10N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    /**
     * @brief This function tests the Triangle2D3N
     * @param rModelPart Model part containing nodes
     * @return true If teh test fails, true otherwise¡
     */
    bool TestTriangle2D3N(ModelPart& rModelPart) {std::stringstream ss; return StreamTestTriangle2D3N(rModelPart, ss);}

    /**
     * @brief This function tests the Triangle2D3N
     * @param rModelPart Model part containing nodes
     * @param rErrorMessage The error message
     * @return true If teh test fails, true otherwise¡
     */
    bool StreamTestTriangle2D3N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    /**
     * @brief This function tests the Triangle2D6N
     * @param rModelPart Model part containing nodes
     * @return true If teh test fails, true otherwise¡
     */
    bool TestTriangle2D6N(ModelPart& rModelPart) {std::stringstream ss; return StreamTestTriangle2D6N(rModelPart, ss);}

    /**
     * @brief This function tests the Triangle2D6N
     * @param rModelPart Model part containing nodes
     * @param rErrorMessage The error message
     * @return true If teh test fails, true otherwise¡
     */
    bool StreamTestTriangle2D6N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    /**
     * @brief This function tests the Quadrilateral2D4N
     * @param rModelPart Model part containing nodes
     * @return true If teh test fails, true otherwise¡
     */
    bool TestQuadrilateral2D4N(ModelPart& rModelPart) {std::stringstream ss; return StreamTestQuadrilateral2D4N(rModelPart, ss);}

    /**
     * @brief This function tests the Quadrilateral2D4N
     * @param rModelPart Model part containing nodes
     * @param rErrorMessage The error message
     * @return true If teh test fails, true otherwise¡
     */
    bool StreamTestQuadrilateral2D4N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    /**
     * @brief This function tests the Quadrilateral2D9N
     * @param rModelPart Model part containing nodes
     * @return true If teh test fails, true otherwise¡
     */
    bool TestQuadrilateral2D9N(ModelPart& rModelPart) {std::stringstream ss; return StreamTestQuadrilateral2D9N(rModelPart, ss);}

    /**
     * @brief This function tests the Quadrilateral2D9N
     * @param rModelPart Model part containing nodes
     * @param rErrorMessage The error message
     * @return true If teh test fails, true otherwise¡
     */
    bool StreamTestQuadrilateral2D9N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    /**
     * @brief This function tests the Quadrilateral2D4N (interface)
     * @param rModelPart Model part containing nodes
     * @return true If teh test fails, true otherwise¡
     */
    bool TestQuadrilateralInterface2D4N(ModelPart& rModelPart) {std::stringstream ss; return StreamTestQuadrilateralInterface2D4N(rModelPart, ss);}

    /**
     * @brief This function tests the Quadrilateral2D4N (interface)
     * @param rModelPart Model part containing nodes
     * @param rErrorMessage The error message
     * @return true If teh test fails, true otherwise¡
     */
    bool StreamTestQuadrilateralInterface2D4N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    /**
     * @brief This function tests the Hexahedra3D8N
     * @param rModelPart Model part containing nodes
     * @return true If teh test fails, true otherwise¡
     */
    bool TestHexahedra3D8N(ModelPart& rModelPart) {std::stringstream ss; return StreamTestHexahedra3D8N(rModelPart, ss);}

    /**
     * @brief This function tests the Hexahedra3D8N
     * @param rModelPart Model part containing nodes
     * @param rErrorMessage The error message
     * @return true If teh test fails, true otherwise¡
     */
    bool StreamTestHexahedra3D8N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    /**
     * @brief This function tests the Hexahedra3D20N
     * @param rModelPart Model part containing nodes
     * @return true If teh test fails, true otherwise¡
     */
    bool TestHexahedra3D20N(ModelPart& rModelPart) {std::stringstream ss; return StreamTestHexahedra3D20N(rModelPart, ss);}

    /**
     * @brief This function tests the Hexahedra3D20N
     * @param rModelPart Model part containing nodes
     * @param rErrorMessage The error message
     * @return true If teh test fails, true otherwise¡
     */
    bool StreamTestHexahedra3D20N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    /**
     * @brief This function tests the Hexahedra3D27N
     * @param rModelPart Model part containing nodes
     * @return true If teh test fails, true otherwise¡
     */
    bool TestHexahedra3D27N(ModelPart& rModelPart) {std::stringstream ss; return StreamTestHexahedra3D27N(rModelPart, ss);}

    /**
     * @brief This function tests the Hexahedra3D27N
     * @param rModelPart Model part containing nodes
     * @param rErrorMessage The error message
     * @return true If teh test fails, true otherwise¡
     */
    bool StreamTestHexahedra3D27N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    /**
     * @brief This function tests the Hexahedra3D8N (interface)
     * @param rModelPart Model part containing nodes
     * @return true If teh test fails, true otherwise¡
     */
    bool TestHexahedraInterface3D8N(ModelPart& rModelPart) {std::stringstream ss; return StreamTestHexahedraInterface3D8N(rModelPart, ss);}

    /**
     * @brief This function tests the Hexahedra3D8N (interface)
     * @param rModelPart Model part containing nodes
     * @param rErrorMessage The error message
     * @return true If teh test fails, true otherwise¡
     */
    bool StreamTestHexahedraInterface3D8N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    /**
     * @brief This function tests the Prism3D6N
     * @param rModelPart Model part containing nodes
     * @return true If teh test fails, true otherwise¡
     */
    bool TestPrism3D6N(ModelPart& rModelPart) {std::stringstream ss; return StreamTestPrism3D6N(rModelPart, ss);}

    /**
     * @brief This function tests the Prism3D6N
     * @param rModelPart Model part containing nodes
     * @param rErrorMessage The error message
     * @return true If teh test fails, true otherwise¡
     */
    bool StreamTestPrism3D6N(
        ModelPart& rModelPart,
        std::stringstream& rErrorMessage
        );

    /**
     * @brief This function tests the Prism3D15N
     * @param rModelPart Model part containing nodes
     * @return true If teh test fails, true otherwise¡
     */
    bool TestPrism3D15N(ModelPart& rModelPart) {std::stringstream ss; return StreamTestPrism3D15N(rModelPart, ss);}

    /**
     * @brief This function tests the Prism3D15N
     * @param rModelPart Model part containing nodes
     * @param rErrorMessage The error message
     * @return true If teh test fails, true otherwise¡
     */
   bool StreamTestPrism3D15N(
       ModelPart& rModelPart,
       std::stringstream& rErrorMessage
       );

    /**
     * @brief This function tests the Prism3D6N (interface)
     * @param rModelPart Model part containing nodes
     * @return true If teh test fails, true otherwise¡
     */
    bool TestPrismInterface3D6N(ModelPart& rModelPart) {std::stringstream ss; return StreamTestPrismInterface3D6N(rModelPart, ss);}

    /**
     * @brief This function tests the Prism3D6N (interface)
     * @param rModelPart Model part containing nodes
     * @param rErrorMessage The error message
     * @return true If teh test fails, true otherwise¡
     */
    bool StreamTestPrismInterface3D6N(
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

    /**
     * @brief This function generates the nodes of the model part
     * @param rModelPart Model part containing nodes
     */
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

    /**
     * @brief Here we verify the area of the element
     * @param rGeometry The geometry
     * @param ThisMethod The integration method
     * @param ReferenceArea The reference area
     * @param rErrorMessage The error message
     */
    bool VerifyAreaByIntegration(
        GeometryType& rGeometry,
        GeometryType::IntegrationMethod ThisMethod,
        const double ReferanceArea,
        std::stringstream& rErrorMessage
        );

    /**
     * @brief Here we verify that a  "displacement field" which varies linearly in space, produces the expected strain distribution.
     * @details This shall be considered a test for shape function derivatives
     * @param rGeometry The geometry
     * @param ThisMethod The integration method
     * @param rErrorMessage The error message
     */
    void VerifyStrainExactness(
        GeometryType& rGeometry,
        GeometryType::IntegrationMethod ThisMethod,
        std::stringstream& rErrorMessage
        );

    /**
     * @brief This function computes the linear strain matrix - useful to verify that a constant strain can be correctly reproduced
     * @param rB The B matrix
     * @param rDN_DX The shape function derivatives
     * @param NumberOfNodes The number of nodes
     * @param Dimension The dimension of the problem
     */
    void CalculateB(
        Matrix& rB,
        Matrix& rDN_DX,
        const std::size_t NumberOfNodes,
        const std::size_t Dimension
        );

    /**
     * @brief This function verifies if the coordinates are inside the geometry
     * @param rGeometry The geometry
     * @param rGlobalCoordinates The global coordinates
     * @param ExpectedResult The expected result
     * @param rErrorMessage The error message
     * @return true if the test passes, false otherwise
     */
    bool VerifyIsInside(
        GeometryType& rGeometry,
        GeometryType::CoordinatesArrayType& rGlobalCoordinates,
        const bool ExpectedResult,
        std::stringstream& rErrorMessage
        );

    /**
     * @brief This function verifies the shape functions
     * @param rGeometry The geometry
     * @param rGlobalCoordinates The global coordinates
     * @param rErrorMessage The error message
     * @return true If teh test fails, true otherwise¡
     */
    bool VerfiyShapeFunctionsValues(
        GeometryType& rGeometry,
        GeometryType::CoordinatesArrayType& rGlobalCoordinates,
        std::stringstream& rErrorMessage
        );

    /**
     * @brief This function verifies the shape functions second derivatives
     * @param rGeometry The geometry
     * @param rGlobalCoordinates The global coordinates
     * @param rErrorMessage The error message
     * @return true If teh test fails, true otherwise¡
     */
    bool VerifyShapeFunctionsSecondDerivativesValues(
        GeometryType& rGeometry,
        GeometryType::CoordinatesArrayType& rGlobalCoordinates,
        std::stringstream& rErrorMessage
        ) const;

    /**
     * @brief This function verifies the shape functions second derivatives interpolations
     * @param rGeometry The geometry
     * @param rErrorMessage The error message
     * @return true If teh test fails, true otherwise¡
     */
    bool VerifyShapeFunctionsSecondDerivativesInterpolation(
    GeometryType& rGeometry,
    std::stringstream& rErrorMessage
    ) const;

    /**
     * @brief Get the name of the intergration method
     * @param rGeometry The geometry
     * @param ThisMethod The integration method
     * @return The integration method name
     */
    std::string GetIntegrationName(
        GeometryType& rGeometry,
        GeometryType::IntegrationMethod ThisMethod
        );

    /**
     * @brief Get the name of the geometry
     * @param rGeometry The geometry
     * @return The geometry name
     */
    std::string GetGeometryName(GeometryType& rGeometry);

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
