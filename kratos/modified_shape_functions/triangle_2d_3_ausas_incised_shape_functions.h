//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Franziska Wahl
//

#if !defined(KRATOS_TRIANGLE_2D_3_AUSAS_INCISED_SHAPE_FUNCTIONS)
#define KRATOS_TRIANGLE_2D_3_AUSAS_INCISED_SHAPE_FUNCTIONS

// System includes

// External includes

// Project includes
#include "modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions.h"

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

class KRATOS_API(KRATOS_CORE) Triangle2D3AusasIncisedShapeFunctions : public Triangle2D3AusasModifiedShapeFunctions
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of Triangle2D3AusasIncisedShapeFunctions
    KRATOS_CLASS_POINTER_DEFINITION(Triangle2D3AusasIncisedShapeFunctions);

    // General type definitions
    typedef AusasModifiedShapeFunctions                        BaseType;
    typedef BaseType::GeometryPointerType                      GeometryPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    Triangle2D3AusasIncisedShapeFunctions(const GeometryPointerType rpInputGeometry,
        const Vector& rNodalDistancesWithExtrapolated, const Vector& rExtrapolatedEdgeRatios);

    /// Destructor
    ~Triangle2D3AusasIncisedShapeFunctions();

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
    ///@name Operations
    ///@{

    /**
    * Returns a reference to the extrapolated edge ratios vector member variable.
    */
    const Vector& GetExtrapolatedEdgeRatios() const ;

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    // Arrays to get edge and node IDs of geometry from edge ID of splitting utility
    const std::array<size_t, 3> edge_id_for_geometry {{2, 0, 1}};
    const std::array<std::array<size_t,2>, 3> node_ids_for_geometry {{{{0,1}}, {{1,2}}, {{2,0}}}};

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void SetPositiveSideCondensationMatrix(Matrix& rPosSideCondMatrix) override;

    void SetNegativeSideCondensationMatrix(Matrix& rNegSideCondMatrix) override;

    /**
    * Returns the intersection points and extrapolated intersection points condensation matrix for
    * positive side Ausas shape functions for incised elements.
    * This matrix is used to transform the subdivisions shape funtion values to the
    * original geometry ones. It has size (n_nodes+n_edges)x(n_nodes).
    * @param rPosSideCondMatrix: Reference to the extrapolated) intersection points condensation matrix to be changed.
    * @param rEdgeNodeI Integers array containing the nodes "I" that conform the edges.
    * @param rEdgeNodeJ Integers array containing the nodes "J" that conform the edges.
    * @param rSplitEdges Integers array containing the original nodes ids and the intersected edges nodes ones.
    */
    void SetPositiveSideCondensationMatrix(
        Matrix& rPosSideCondMatrix,
        const std::vector<int>& rEdgeNodeI,
        const std::vector<int>& rEdgeNodeJ,
        const std::vector<int>& rSplitEdges);

    /**
    * Returns the intersection points and extrapolated intersection points condensation matrix for
    * negative side Ausas shape functions for incised elements.
    * This matrix is used to transform the subdivisions shape funtion values to the
    * original geometry ones. It has size (n_nodes+n_edges)x(n_nodes).
    * @param rNegSideCondMatrix: Reference to the (extrapolated) intersection points condensation matrix to be changed.
    * @param rEdgeNodeI Integers array containing the nodes "I" that conform the edges.
    * @param rEdgeNodeJ Integers array containing the nodes "J" that conform the edges.
    * @param rSplitEdges Integers array containing the original nodes ids and the intersected edges nodes ones.
    */
    void SetNegativeSideCondensationMatrix(
        Matrix& rNegSideCondMatrix,
        const std::vector<int>& rEdgeNodeI,
        const std::vector<int>& rEdgeNodeJ,
        const std::vector<int>& rSplitEdges);

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

    const Vector mExtraEdgeRatios;

    ///@}
    ///@name Serialization
    ///@{

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

    /// Assignment operator.
    Triangle2D3AusasIncisedShapeFunctions& operator=(Triangle2D3AusasModifiedShapeFunctions const& rOther);

    /// Copy constructor.
    Triangle2D3AusasIncisedShapeFunctions(Triangle2D3AusasIncisedShapeFunctions const& rOther) :
        Triangle2D3AusasModifiedShapeFunctions(rOther.GetInputGeometry(), rOther.GetNodalDistances()),
        mExtraEdgeRatios(rOther.mExtraEdgeRatios) {
    };

    ///@}

};// class Triangle2D3AusasIncisedShapeFunctions

}
#endif /* KRATOS_TRIANGLE_2D_3_AUSAS_INCISED_SHAPE_FUNCTIONS defined */
