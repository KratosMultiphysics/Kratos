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

    /**
    * Returns the shape function values in the positive split element side for a given quadrature.
    * @return rPositiveSideShapeFunctionValues: Matrix containing the positive side computed shape function values.
    * @return rPositiveSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the positive side.
    * @return rPositiveSideWeightsValues: Vector containing the Gauss pts. positive side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputePositiveSideShapeFunctionsAndGradientsValues(
        Matrix &rPositiveSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rPositiveSideShapeFunctionsGradientsValues,
        Vector &rPositiveSideWeightsValues,
        const IntegrationMethodType IntegrationMethod) override;

    /**
    * Returns the shape function values in the negative split element side for a given quadrature.
    * @return rNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeNegativeSideShapeFunctionsAndGradientsValues(
        Matrix &rNegativeSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rNegativeSideShapeFunctionsGradientsValues,
        Vector &rNegativeSideWeightsValues,
        const IntegrationMethodType IntegrationMethod) override;

    /**
    * Returns the shape function values in the positive split element interface side for a given quadrature.
    * @return rInterfacePositiveSideShapeFunctionValues: Matrix containing the positive side computed shape function values.
    * @return rInterfacePositiveSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the positive side.
    * @return rInterfacePositiveSideWeightsValues: Vector containing the Gauss pts. positive side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
        Matrix &rInterfacePositiveSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rInterfacePositiveSideShapeFunctionsGradientsValues,
        Vector &rInterfacePositiveSideWeightsValues,
        const IntegrationMethodType IntegrationMethod) override;

    /**
    * Returns the shape function values in the negative split element side for a given quadrature.
    * @return rInterfaceNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rInterfaceNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rInterfaceNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
        Matrix &rInterfaceNegativeSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rInterfaceNegativeSideShapeFunctionsGradientsValues,
        Vector &rInterfaceNegativeSideWeightsValues,
        const IntegrationMethodType IntegrationMethod) override;

    /**
    * Given a face id, returns the shape function values in the positive split element exterior face side for a given quadrature.
    * @return rInterfacePositiveSideShapeFunctionValues: Matrix containing the positive side computed shape function values.
    * @return rInterfacePositiveSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the positive side.
    * @return rInterfacePositiveSideWeightsValues: Vector containing the Gauss pts. positive side weights (already multiplied by the Jacobian).
    * @param FaceId Face local id. in where the values are to be computed.
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
        Matrix &rPositiveExteriorFaceShapeFunctionsValues,
        ShapeFunctionsGradientsType &rPositiveExteriorFaceShapeFunctionsGradientsValues,
        Vector &rPositiveExteriorFaceWeightsValues,
        const unsigned int FaceId,
        const IntegrationMethodType IntegrationMethod) override;

    /**
    * Given a face id, returns the shape function values in the negative split element exterior face side for a given quadrature.
    * @return rInterfaceNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rInterfaceNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rInterfaceNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param FaceId Face local id. in where the values are to be computed.
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
        Matrix &rNegativeExteriorFaceShapeFunctionsValues,
        ShapeFunctionsGradientsType &rNegativeExteriorFaceShapeFunctionsGradientsValues,
        Vector &rNegativeExteriorFaceWeightsValues,
        const unsigned int FaceId,
        const IntegrationMethodType IntegrationMethod) override;

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

    /**
     * @brief Set the Positive Side Condensation Matrix object
     * This function sets the positive side condensation matrix
     * Note that for the Ausas FE space a different condensation matrix is required for the positive and negative sides
     * @param rPosSideCondMatrix Positive side condensation matrix
     */
    virtual void SetPositiveSideCondensationMatrix(Matrix& rPosSideCondMatrix) = 0;

    /**
     * @brief Set the Negative Side Condensation Matrix object
     * This function sets the negative side condensation matrix
     * Note that for the Ausas FE space a different condensation matrix is required for the positive and negative sides
     * @param rNegSideCondMatrix Negative side condensation matrix
     */
    virtual void SetNegativeSideCondensationMatrix(Matrix& rNegSideCondMatrix) = 0;

    /**
    * Returns the intersection points condensation matrix for positive side Ausas sh functions.
    * This matrix is used to extrapolate the subdivisions shape funtion values to the
    * original geometry ones. It has size (nnodes+nedges)x(nnodes).
    * @return rPosSideCondMatrix: Reference to the intersection points condensation matrix.
    * @param rEdgeNodeI Integers array containing the nodes "I" that conform the edges.
    * @param rEdgeNodeJ Integers array containing the nodes "J" that conform the edges.
    * @param rSplitEdges Integers array containing the original nodes ids and the intersected edges nodes ones.
    */
    virtual void SetPositiveSideCondensationMatrix(
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
    virtual void SetNegativeSideCondensationMatrix(
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
