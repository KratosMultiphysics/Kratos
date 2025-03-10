//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_AUSAS_MODIFIED_SHAPE_FUNCTIONS)
#define KRATOS_AUSAS_MODIFIED_SHAPE_FUNCTIONS

// System includes

// External includes

// Project includes
#include "modified_shape_functions/modified_shape_functions.h"

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

class KRATOS_API(KRATOS_CORE) AusasModifiedShapeFunctions : public ModifiedShapeFunctions
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of AusasModifiedShapeFunctions
    KRATOS_CLASS_POINTER_DEFINITION(AusasModifiedShapeFunctions);

    // General type definitions
    typedef ModifiedShapeFunctions                             BaseType;
    typedef BaseType::GeometryType                             GeometryType;
    typedef BaseType::GeometryPointerType                      GeometryPointerType;
    typedef BaseType::IntegrationMethodType                    IntegrationMethodType;
    typedef BaseType::ShapeFunctionsGradientsType              ShapeFunctionsGradientsType;

    typedef BaseType::IndexedPointGeometryType                 IndexedPointGeometryType;
    typedef BaseType::IndexedPointGeometryPointerType          IndexedPointGeometryPointerType;

    typedef BaseType::IntegrationPointType                     IntegrationPointType;
    typedef BaseType::IntegrationPointsArrayType               IntegrationPointsArrayType;
    typedef BaseType::IntegrationPointsContainerType           IntegrationPointsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    AusasModifiedShapeFunctions(const GeometryPointerType rpInputGeometry, const Vector& rNodalDistances);

    /// Destructor
    ~AusasModifiedShapeFunctions();

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

    void SetCondensationMatrix(Matrix& rIntPointCondMatrix) override
    {
        KRATOS_ERROR << "\'SetCondensationMatrix\' cannot be used with the Ausas FE space. Use either \'SetPositiveSideCondensationMatrix\' or \'SetNegativeSideCondensationMatrix\' instead." << std::endl;
    }

    void SetPositiveSideCondensationMatrix(Matrix& rPosSideCondMatrix) override
    {
        KRATOS_ERROR << "Calling Ausas base class \'SetPositiveSideCondensationMatrix\'. Call the derived class one instead." << std::endl;
    }

    void SetNegativeSideCondensationMatrix(Matrix& rNegSideCondMatrix) override
    {
        KRATOS_ERROR << "Calling Ausas base class \'SetNegativeSideCondensationMatrix\'. Call the derived class one instead." << std::endl;
    }

    /**
    * Returns the intersection points condensation matrix for positive side Ausas sh functions.
    * This matrix is used to extrapolate the subdivisions shape funtion values to the
    * original geometry ones. It has size (nnodes+nedges)x(nnodes).
    * @return rPosSideCondMatrix: Reference to the intersection points condensation matrix.
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
    * Returns the intersection points condensation matrix for negative side Ausas sh functions.
    * This matrix is used to extrapolate the subdivisions shape funtion values to the
    * original geometry ones. It has size (nnodes+nedges)x(nnodes).
    * @return rNegSideCondMatrix: Reference to the intersection points condensation matrix.
    * @param rEdgeNodeI Integers array containing the nodes "I" that conform the edges.
    * @param rEdgeNodeJ Integers array containing the nodes "J" that conform the edges.
    * @param rSplitEdges Integers array containing the original nodes ids and the intersected edges nodes ones.
    */
    void SetNegativeSideCondensationMatrix(
        Matrix& rNegSideCondMatrix,
        const std::vector<int>& rEdgeNodeI,
        const std::vector<int>& rEdgeNodeJ,
        const std::vector<int>& rSplitEdges);
       
    /**
    * Returns the negative side tangential vector (in the direction of contact line).
    * @param FaceIndices: vector giving the indices of the contact faces from DivideGeometry::mContactFace
    * @return rNegativeSideContactLineVector: single vector showing the contact line. 
    */
    void ComputeNegativeSideContactLineVector(
        std::vector<unsigned int>& FaceIndices,
        std::vector<Vector> &rNegativeSideContactLineVector) override
    {
        FaceIndices.clear();
        rNegativeSideContactLineVector.clear();
    }
    /* bool ComputeNegativeSideContactLineVector(
        Vector &rNegativeSideContactLineVector) override
    {
        rNegativeSideContactLineVector = ZeroVector(3);
        return false;
    } */

    /**
    * Returns the shape function values in the negative split element side for a given quadrature on the contact line.
    * @param ContactLineIndices: indices associated with the outer faces that can be considered as contact lines
    * @return rContactLineNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rContactLineNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rContactLineNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeContactLineNegativeSideShapeFunctionsAndGradientsValues(
        std::vector<unsigned int>& ContactLineIndices,
        std::vector<Matrix> &rContactLineNegativeSideShapeFunctionsValues,
        std::vector<ShapeFunctionsGradientsType> &rContactLineNegativeSideShapeFunctionsGradientsValues,
        std::vector<Vector> &rContactLineNegativeSideWeightsValues,
        const IntegrationMethodType IntegrationMethod) override
    {
        double dummy = 0.0;
    }
    /* void ComputeContactLineNegativeSideShapeFunctionsAndGradientsValues(
        Matrix &rContactLineNegativeSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rContactLineNegativeSideShapeFunctionsGradientsValues,
        Vector &rContactLineNegativeSideWeightsValues,
        const IntegrationMethodType IntegrationMethod) override */ 

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
    AusasModifiedShapeFunctions& operator=(AusasModifiedShapeFunctions const& rOther);

    /// Copy constructor.
    AusasModifiedShapeFunctions(AusasModifiedShapeFunctions const& rOther) :
        ModifiedShapeFunctions(rOther.GetInputGeometry(), rOther.GetNodalDistances()) {
    };

    ///@}

};// class AusasModifiedShapeFunctions

}
#endif /* KRATOS_AUSAS_MODIFIED_SHAPE_FUNCTIONS defined */
