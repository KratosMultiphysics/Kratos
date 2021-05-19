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

#if !defined(KRATOS_TETRAHEDRA_3D_4_AUSAS_INCISED_SHAPE_FUNCTIONS)
#define KRATOS_TETRAHEDRA_3D_4_AUSAS_INCISED_SHAPE_FUNCTIONS

// System includes

// External includes

// Project includes
#include "modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions.h"

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

class KRATOS_API(KRATOS_CORE) Tetrahedra3D4AusasIncisedShapeFunctions : public Tetrahedra3D4AusasModifiedShapeFunctions
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of Tetrahedra3D4AusasModifiedShapeFunctions
    KRATOS_CLASS_POINTER_DEFINITION(Tetrahedra3D4AusasIncisedShapeFunctions);

    // General type definitions
    typedef AusasModifiedShapeFunctions                        BaseType;
    typedef BaseType::GeometryPointerType                      GeometryPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    Tetrahedra3D4AusasIncisedShapeFunctions(
        const GeometryPointerType rpInputGeometry,
        const Vector& rNodalDistancesWithExtrapolated,
        const Vector& rExtrapolatedEdgeRatios);

    /// Destructor
    ~Tetrahedra3D4AusasIncisedShapeFunctions();

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
    const std::array<size_t, 6> edge_id_for_geometry {{0, 2, 3, 1, 4, 5}};
    const std::array<std::array<size_t,2>, 6> node_ids_for_geometry {{{{0,1}}, {{2,0}}, {{0,3}}, {{1,2}}, {{1,3}}, {{2,3}}}};

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
    Tetrahedra3D4AusasIncisedShapeFunctions& operator=(Tetrahedra3D4AusasIncisedShapeFunctions const& rOther);

    /// Copy constructor.
    Tetrahedra3D4AusasIncisedShapeFunctions(Tetrahedra3D4AusasIncisedShapeFunctions const& rOther) :
        Tetrahedra3D4AusasModifiedShapeFunctions(rOther.GetInputGeometry(), rOther.GetNodalDistances()),
        mExtraEdgeRatios(rOther.mExtraEdgeRatios) {
    };

    ///@}

};// class Tetrahedra3D4AusasIncisedShapeFunctions

}
#endif /* KRATOS_TETRAHEDRA_3D_4_AUSAS_INCISED_SHAPE_FUNCTIONS defined */
