// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

#if !defined(KRATOS_SBM_LAPLACIAN_CONDITION_H_INCLUDED)
#define  KRATOS_SBM_LAPLACIAN_CONDITION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "includes/condition.h"
#include "geometries/geometry.h"
#include "includes/variables.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) SBMLaplacianCondition: public Condition
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of SBMLaplacianCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SBMLaplacianCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    SBMLaplacianCondition(
        IndexType NewId,
        Geometry< Node >::Pointer pGeometry);

    SBMLaplacianCondition(
        IndexType NewId,
        Geometry< Node >::Pointer pGeometry,
        Properties::Pointer pProperties);

    SBMLaplacianCondition(
        IndexType NewId,
        const NodesArrayType& ThisNodes);

    /// Destructor.
    ~SBMLaplacianCondition() override;

    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        Properties::Pointer pProperties) const override;

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        Properties::Pointer pProperties) const override;

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(
        DofsVectorType& ConditionalDofList,
        const ProcessInfo& CurrentProcessInfo) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

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

protected:

    ///@name Protected Life Cycle
    ///@{

    // Internal default constructor for serialization
    SBMLaplacianCondition();

    ///@}
    ///@name Protected Operations
    ///@{


    ///@}

private:

    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    SBMLaplacianCondition& operator=(SBMLaplacianCondition const& rOther);

    /// Copy constructor.
    SBMLaplacianCondition(SBMLaplacianCondition const& rOther);

    ///@}

}; // Class SBMLaplacianCondition

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_LAPLACIAN_SHIFTED_BOUNDARY_CONDITION_H_INCLUDED  defined