//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

#ifndef KRATOS_MPI_NORMAL_CALCULATION_UTILITIES_H
#define KRATOS_MPI_NORMAL_CALCULATION_UTILITIES_H

// System includes

// External includes


// Project includes
#include "includes/model_part.h"


namespace Kratos
{
///@addtogroup KratosCore
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

/// Some tools to calculate face and nodal normals on an MPI partitioned environment
class KRATOS_API(KRATOS_MPI_CORE) MPINormalCalculationUtils
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MPINormalCalculationUtils
    KRATOS_CLASS_POINTER_DEFINITION(MPINormalCalculationUtils);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MPINormalCalculationUtils();

    /// Destructor.
    virtual ~MPINormalCalculationUtils();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    int Check(ModelPart& rModelPart);

    void OrientFaces(ModelPart& rModelPart,
                     bool OutwardsPositive);

    void CalculateOnSimplex(ModelPart& rModelPart,
                            int Dimension,
                            const Variable<double>& rVariable,
                            const double rAlpha);

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
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;


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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    void IdentifyFaces(ModelPart& rModelPart,
                       const Variable<double>& rVariable,
                       int& MaxNeigh,
                       int& NumNodes);


    void InitializeNormalData(ModelPart& rModelPart,
                              const Variable<double>& rVariable,
                              int MaxNeigh,
                              std::vector<double>& rNormals,
                              std::vector<int>& pActiveNeigh,
                              const int NumNodes);

    void DetectEdges(ModelPart& rModelPart,
                     unsigned int Dimension,
                     double MaxAngle,
                     const std::vector<double>& rNormals,
                     const std::vector<int>& rActiveNeigh,
                     const int NumNodes,
                     const int MaxNeighs);


    bool OrientElement(Geometry< Node >& rGeom);


    void NormalContribution(Geometry< Node >&rGeom);


    void FaceNormal2D(array_1d<double,3>& An,
                      Geometry<Node >& rGeometry);


    void FaceNormal3D(array_1d<double,3>& An,
                      Geometry<Node >& rGeometry);

    void UpdateNodeNormals(ModelPart& rModelPart,
                           const unsigned int Dimension,
                           const Variable<double>& rVariable);


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
    MPINormalCalculationUtils& operator=(MPINormalCalculationUtils const& rOther);

    /// Copy constructor.
    MPINormalCalculationUtils(MPINormalCalculationUtils const& rOther);


    ///@}

}; // Class MPINormalCalculationUtils

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MPINormalCalculationUtils& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MPINormalCalculationUtils& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MPI_NORMAL_CALCULATION_UTILITIES_H
