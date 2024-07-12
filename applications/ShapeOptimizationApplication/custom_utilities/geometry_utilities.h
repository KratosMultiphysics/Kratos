// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef GEOMETRY_UTILITIES_H
#define GEOMETRY_UTILITIES_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// ==============================================================================

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

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.

*/

class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) GeometryUtilities
{
public:
    ///@name Type Definitions
    ///@{

    // For better reading
    typedef array_1d<double,3> array_3d;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    typedef ModelPart::ElementType::GeometryType GeometryType;
    typedef std::size_t SizeType;
    using NodeType = Node;

    /// Pointer definition of GeometryUtilities
    KRATOS_CLASS_POINTER_DEFINITION(GeometryUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GeometryUtilities( ModelPart& modelPart )
        : mrModelPart( modelPart )
    {
    }

    /// Destructor.
    virtual ~GeometryUtilities()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void CalculateNodalAreasFromConditions();

    void ComputeUnitSurfaceNormals();

    void ProjectNodalVariableOnUnitSurfaceNormals( const Variable<array_3d> &rNodalVariable );

    void ProjectNodalVariableOnDirection( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rDirectionVariable);

    void ProjectNodalVariableOnTangentPlane( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rPlaneNormalVariable);

    void ExtractBoundaryNodes( std::string const& rBoundarySubModelPartName );

    void ExtractEdgeNodes( std::string const& rEdgeSubModelPartName );

    std::tuple<std::vector<double>, std::vector<double>> ComputeDistancesToBoundingModelPart(ModelPart& rBoundingModelPart);

    template <class TContainerType>
    double CalculateLength(TContainerType& rContainer)
    {
        double length = block_for_each<SumReduction<double>>(rContainer, [&](typename TContainerType::value_type& rEntity){
            return rEntity.GetGeometry().Length();
        });

        return length;
    }

    double ComputeVolume();

    void ComputeVolumeShapeDerivatives(const Variable<array_3d>& rDerivativeVariable);

    // Computes the gaussian curvature at all surface nodes depending on the surrounding elements
    // (I) Quadratic Elements:  Curvature tensor by shape functions
    // (II) Linear Quadrangle:  Taubin http://graphics.stanford.edu/courses/cs348a-17-winter/ReaderNotes/taubin-iccv95b.pdf
    //                          Camprubi Estebo https://mediatum.ub.tum.de/doc/601055/00000039.pdf
    // (III) Linear Triangle:   Meyer https://authors.library.caltech.edu/99186/2/diffGeoOps.pdf
    void CalculateGaussianCurvature();

    // Computes the gaussian curvature at a surface node from the curvature tensor.
    // For quadratic elements (e.g. 6NTriangle).
    double GaussianCurvatureForNodeFromTensor(const NodeType &rNode);

    // Computes the gaussian curvature at a surface node according to Taubin.
    // For 4NQuads, variable naming according to Camprubi Estebo.
    double GaussianCurvatureForNodeTaubin(const NodeType &rNode);

    // Computes the gaussian curvature at a surface node according to Meyer.
    // For 3NTriangles.
    double GaussianCurvatureForNodeMeyer(const NodeType &rNode);

    // Selects the curvature technique according to the surrounding elements
    // Precedence: (I) > (II) > (III)
    std::string GetCurvatureTechnique(const NodeType &rNode);

    bool CheckIfElementIsQuadratic(const Kratos::GlobalPointer<Kratos::Condition> pElement);

    bool CheckIfNodesHasQuadraticNeigbourElement(const NodeType &rNode);

    void LocalPointInElement(const NodeType& rNode, const Kratos::GlobalPointer<Kratos::Condition> pElement, Kratos::Point::CoordinatesArrayType& rLocalPoint);

    Matrix CurvatureTensor(const NodeType& rNode, const Kratos::GlobalPointer<Kratos::Condition> pElement);

    void BaseVectors(const NodeType& rNode, const Kratos::GlobalPointer<Kratos::Condition> pElement, Vector& rG1, Vector& rG2);

    void CartesianBaseVectors(const NodeType& rNode, const Kratos::GlobalPointer<Kratos::Condition> pElement, Vector& rE1, Vector& rE2);

    void TransformTensorCoefficients(Matrix& rTensor, Matrix& rResultTensor, Vector&rG1, Vector&rG2, Vector&rE1, Vector&rE2);

    void InnerAngleAndMixedAreaOf3D3NTriangletAtNode(const NodeType& rNode, const Kratos::GlobalPointer<Kratos::Condition> pElement, double& rInnerAngle, double& rMixedArea);

    // --------------------------------------------------------------------------

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
        return "GeometryUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GeometryUtilities";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


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

    ModelPart& mrModelPart;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // --------------------------------------------------------------------------
    void CalculateAreaNormalsFromConditions();

    // --------------------------------------------------------------------------
    void CalculateUnitNormals();

    // --------------------------------------------------------------------------


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
//      GeometryUtilities& operator=(GeometryUtilities const& rOther);

    /// Copy constructor.
//      GeometryUtilities(GeometryUtilities const& rOther);


    ///@}

}; // Class GeometryUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // GEOMETRY_UTILITIES_H
