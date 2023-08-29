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

#if !defined(KRATOS_TETRAHEDRA_3D_4_MODIFIED_SHAPE_FUNCTIONS)
#define KRATOS_TETRAHEDRA_3D_4_MODIFIED_SHAPE_FUNCTIONS

// System includes

// External includes

// Project includes
#include "utilities/divide_tetrahedra_3d_4.h"
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

class KRATOS_API(KRATOS_CORE) Tetrahedra3D4ModifiedShapeFunctions : public ModifiedShapeFunctions
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of Tetrahedra3D4ModifiedShapeFunctions
    KRATOS_CLASS_POINTER_DEFINITION(Tetrahedra3D4ModifiedShapeFunctions);

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
    Tetrahedra3D4ModifiedShapeFunctions(const GeometryPointerType rpInputGeometry, const Vector& rNodalDistances);

    /// Destructor
    ~Tetrahedra3D4ModifiedShapeFunctions() override;

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
    * Returns the member pointer to the splitting utility.
    */
    const DivideGeometry<Node>::Pointer pGetSplittingUtil() const override;

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

    void SetCondensationMatrix(Matrix& rIntPointCondMatrix) override;

    void SetPositiveSideCondensationMatrix(Matrix& rPosSideCondMatrix) override;

    void SetNegativeSideCondensationMatrix(Matrix& rNegSideCondMatrix) override;

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

    DivideTetrahedra3D4<Node>::Pointer mpTetrahedraSplitter;

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
    Tetrahedra3D4ModifiedShapeFunctions& operator=(Tetrahedra3D4ModifiedShapeFunctions const& rOther);

    /// Copy constructor.
    Tetrahedra3D4ModifiedShapeFunctions(Tetrahedra3D4ModifiedShapeFunctions const& rOther) :
        ModifiedShapeFunctions(rOther.GetInputGeometry(), rOther.GetNodalDistances()),
        mpTetrahedraSplitter(new DivideTetrahedra3D4<Node>(*rOther.GetInputGeometry(), rOther.GetNodalDistances())) {

        // Perform the element splitting
        mpTetrahedraSplitter->GenerateDivision();
        mpTetrahedraSplitter->GenerateIntersectionsSkin();
    };

    ///@}

};// class Tetrahedra3D4ModifiedShapeFunctions

}
#endif /* KRATOS_TETRAHEDRA_3D_4_MODIFIED_SHAPE_FUNCTIONS defined */
