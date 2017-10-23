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

#if !defined(KRATOS_TRIANGLE_2D_3_MODIFIED_SHAPE_FUNCTIONS)
#define KRATOS_TRIANGLE_2D_3_MODIFIED_SHAPE_FUNCTIONS

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

class KRATOS_API(KRATOS_CORE) Triangle2D3ModifiedShapeFunctions : public ModifiedShapeFunctions
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of Triangle2D3ModifiedShapeFunctions
    KRATOS_CLASS_POINTER_DEFINITION(Triangle2D3ModifiedShapeFunctions);

    // General type definitions
    typedef ModifiedShapeFunctions                                      BaseType;
    typedef typename BaseType::GeometryType                             GeometryType;
    typedef typename BaseType::GeometryPointerType                             GeometryPointerType;
    typedef typename BaseType::IntegrationMethodType                    IntegrationMethodType;
    typedef typename BaseType::ShapeFunctionsGradientsType              ShapeFunctionsGradientsType;
    
    typedef typename BaseType::IndexedPointGeometryType                 IndexedPointGeometryType;
    typedef typename BaseType::IndexedPointGeometryPointerType          IndexedPointGeometryPointerType;

    typedef typename BaseType::IntegrationPointType                     IntegrationPointType;
    typedef typename BaseType::IntegrationPointsArrayType               IntegrationPointsArrayType;
    typedef typename BaseType::IntegrationPointsContainerType           IntegrationPointsContainerType;

    int mSplitEdgesNumber;  // Number of split edges
    int mDivisionsNumber;   // Number of generated subdivisions

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    Triangle2D3ModifiedShapeFunctions(GeometryPointerType rpInputGeometry, Vector& rNodalDistances);

    /// Destructor
    ~Triangle2D3ModifiedShapeFunctions();

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
    * Returns the shape function values in both the positive or negative split element sides for a given quadrature.
    * @return rPositiveSideShapeFunctionValues: Matrix containing the positive side computed shape function values.
    * @return rNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rPositiveSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the positive side.
    * @return rNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rPositiveSideWeightsValues: Vector containing the Gauss pts. positive side weights (already multiplied by the Jacobian).
    * @return rNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod: Desired integration quadrature.
    */
    void GetShapeFunctionsAndGradientsValues(Matrix &rPositiveSideShapeFunctionsValues,
                                             Matrix &rNegativeSideShapeFunctionsValues,
                                             std::vector<Matrix> &rPositiveSideShapeFunctionsGradientsValues,
                                             std::vector<Matrix> &rNegativeSideShapeFunctionsGradientsValues,
                                             Vector &rPositiveSideWeightsValues,
                                             Vector &rNegativeSideWeightsValues,
                                             const IntegrationMethodType IntegrationMethod) override;

    /**
    * Returns the shape function values in the positive split element side for a given quadrature.
    * @return rPositiveSideShapeFunctionValues: Matrix containing the positive side computed shape function values.
    * @return rPositiveSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the positive side.
    * @return rPositiveSideWeightsValues: Vector containing the Gauss pts. positive side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod: Desired integration quadrature.
    */
    void GetPositiveSideShapeFunctionsAndGradientsValues(Matrix &rPositiveSideShapeFunctionsValues,
                                                         std::vector<Matrix> &rPositiveSideShapeFunctionsGradientsValues,
                                                         Vector &rPositiveSideWeightsValues,
                                                         const IntegrationMethodType IntegrationMethod) override;

    /**
    * Returns the shape function values in the negative split element side for a given quadrature.
    * @return rNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod: Desired integration quadrature.
    */
    void GetNegativeSideShapeFunctionsAndGradientsValues(Matrix &rNegativeSideShapeFunctionsValues,
                                                         std::vector<Matrix> &rNegativeSideShapeFunctionsGradientsValues,
                                                         Vector &rNegativeSideWeightsValues,
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

    /**
    * Returns the shape function values in either the positive or negative element subdivision for a given quadrature.
    * @return rShapeFunctionValues: Matrix containing the computed shape function values.
    * @return rShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values.
    * @return rWeightsValues: Vector containing the Gauss pts. weights (already multiplied by the Jacobian).
    * @param rSubdivisionGeom: std::vector of subdivisions point based geometries where the values are to be computed.
    * @param IntegrationMethod: Desired integration quadrature.
    */
    void ComputeValuesOnOneSide(Matrix &rShapeFunctionsValues,
                                std::vector<Matrix> &rShapeFunctionsGradientsValues,
                                Vector &rWeightsValues,
                                const std::vector<IndexedPointGeometryPointerType> &rSubdivisionsVector,
                                const Matrix &p_matrix,
                                const IntegrationMethodType IntegrationMethod) override;

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

    // /// Assignment operator.
    // Triangle2D3ModifiedShapeFunctions& operator=(Triangle2D3ModifiedShapeFunctions const& rOther);

    // /// Copy constructor.
    // Triangle2D3ModifiedShapeFunctions(Triangle2D3ModifiedShapeFunctions const& rOther)
    //     : mrInputGeometry(rOther.mrInputGeometry) , mrNodalDistances(rOther.mrNodalDistances) {};

    ///@}

};// class Triangle2D3ModifiedShapeFunctions

}
#endif /* KRATOS_TRIANGLE_2D_3_MODIFIED_SHAPE_FUNCTIONS defined */
