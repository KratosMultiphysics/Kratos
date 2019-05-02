//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Guillermo Casas gcasas@gmail.com
//

#ifndef KRATOS_MONOLITHIC_DEM_COUPLED_WALL_CONDITION_H
#define KRATOS_MONOLITHIC_DEM_COUPLED_WALL_CONDITION_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/process_info.h"

// Application includes
#include "../FluidDynamicsApplication/fluid_dynamics_application_variables.h"
#include "includes/deprecated_variables.h"
#include "includes/cfd_variables.h"

//Other Applications includes
#include "../../applications/FluidDynamicsApplication/fluid_dynamics_application.h"
#include "../../applications/FluidDynamicsApplication/custom_conditions/monolithic_wall_condition.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
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

/// Implements a wall condition for the monolithic formulation.
/**
  It is intended to be used in combination with ASGS and VMS elements or their derived classes
  and the ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent time scheme, which supports
  slip conditions.
  This condition will add a wall stress term to all nodes identified with IS_STRUCTURE!=0.0 (in the
  non-historic database, that is, assigned using Node.SetValue()). This stress term is determined
  according to the wall distance provided as Y_WALL.
  @see ASGS2D,ASGS3D,VMS,ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent
 */
template< unsigned int TDim, unsigned int TNumNodes = TDim >
class KRATOS_API(SWIMMING_DEM_APPLICATION) MonolithicDEMCoupledWallCondition : public MonolithicWallCondition < TDim, TNumNodes >
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(MonolithicDEMCoupledWallCondition);

    typedef Node < 3 > NodeType;

    typedef Properties PropertiesType;

    typedef Geometry<NodeType> GeometryType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    /// Pointer definition of MonolithicDEMCoupledWallCondition

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new condition
      */
    MonolithicDEMCoupledWallCondition(IndexType NewId = 0):
        MonolithicWallCondition< TDim, TNumNodes >(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    MonolithicDEMCoupledWallCondition(IndexType NewId,
                            const NodesArrayType& ThisNodes):
        MonolithicWallCondition< TDim, TNumNodes >(NewId,ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    MonolithicDEMCoupledWallCondition(IndexType NewId,
                            GeometryType::Pointer pGeometry):
        MonolithicWallCondition< TDim, TNumNodes >(NewId,pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    MonolithicDEMCoupledWallCondition(IndexType NewId,
                            GeometryType::Pointer pGeometry,
                            PropertiesType::Pointer pProperties):
        MonolithicWallCondition< TDim, TNumNodes >(NewId,pGeometry,pProperties)
    {
    }

    /// Copy constructor.
    MonolithicDEMCoupledWallCondition(MonolithicDEMCoupledWallCondition const& rOther):
        MonolithicWallCondition< TDim, TNumNodes >(rOther)
    {
    }

    /// Destructor.
    virtual ~MonolithicDEMCoupledWallCondition() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    MonolithicDEMCoupledWallCondition & operator=(MonolithicDEMCoupledWallCondition const& rOther)
    {
        Condition::operator=(rOther);

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new MonolithicDEMCoupledWallCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Condition::Pointer(new MonolithicDEMCoupledWallCondition(NewId, Condition::GetGeometry().Create(ThisNodes), pProperties));
    }
    Condition::Pointer Create(IndexType NewId,
                           GeometryType::Pointer pGeom,
                           PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared< MonolithicDEMCoupledWallCondition >(NewId, pGeom, pProperties);
    }

    /// Return local contributions of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see CalculateLocalVelocityContribution
      */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType BlockSize = (rCurrentProcessInfo[FRACTIONAL_STEP] == 1) ? TDim + 1 : TDim;
        const SizeType LocalSize = BlockSize * TNumNodes;

        if (rLeftHandSideMatrix.size1() != LocalSize)
            rLeftHandSideMatrix.resize(LocalSize,LocalSize);

        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);
        noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    }

    /// Return a matrix of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see DampingMatrix
      */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType BlockSize = (rCurrentProcessInfo[FRACTIONAL_STEP] == 1) ? TDim + 1 : TDim;
        const SizeType LocalSize = BlockSize * TNumNodes;

        if (rLeftHandSideMatrix.size1() != LocalSize)
            rLeftHandSideMatrix.resize(LocalSize,LocalSize);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);
    }

    /// Return local right hand side of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see CalculateLocalVelocityContribution
      */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo) override
    {
        if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1){
            const SizeType BlockSize = TDim + 1;
            const SizeType LocalSize = BlockSize * TNumNodes;

            if (rRightHandSideVector.size() != LocalSize)
                rRightHandSideVector.resize(LocalSize);

            noalias(rRightHandSideVector) = ZeroVector(LocalSize);
        }

        else {
            const unsigned int BlockSize = TDim;
            const SizeType LocalSize = BlockSize * TNumNodes;

            if (rRightHandSideVector.size() != LocalSize)
                rRightHandSideVector.resize(LocalSize);

            noalias(rRightHandSideVector) = ZeroVector(LocalSize);

//            const unsigned int LocalSize = TDim;

//            const GeometryType& rGeom = this->GetGeometry();
//            const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
//            const unsigned int NumGauss = IntegrationPoints.size();

//            MatrixType NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

//            array_1d<double,3> Normal;
//            this->CalculateNormal(Normal); //this already contains the area
//            double A = std::sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]+Normal[2]*Normal[2]);
//            Normal /= A;

//            // CAUTION: "Jacobian" is 2.0*A for triangles but 0.5*A for lines
//            double J = (TDim == 2) ? 0.5*A : 2.0*A;

//            for (unsigned int g = 0; g < NumGauss; g++)
//            {
//                Vector N = row(NContainer,g);
//                double Weight = J * IntegrationPoints[g].Weight();

//                // Neumann boundary condition
//                for (unsigned int i = 0; i < TNumNodes; i++)
//                {
//                    //unsigned int row = i*LocalSize;
//                    const NodeType& rConstNode = this->GetGeometry()[i];
//                    const array_1d<double,3>& rVel = this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY);

//                    for (unsigned int j = 0; j < TNumNodes; j++)
//                    {
//                        unsigned int row = j*LocalSize;
//                        for (unsigned int d = 0; d < TDim;d++)
//                            rLocalVector[row+d] -= Weight*DN[j]*N[i]*rVel[d]*Normal[d];
//                    }
//                }
//            }
        }
    }


    /// Provides the global indices for each one of this element's local rows.
    /** This determines the elemental equation ID vector for all elemental DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    void EquationIdVector(EquationIdVectorType& rResult,
                                  ProcessInfo& rCurrentProcessInfo) override;


    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& ConditionDofList,
                            ProcessInfo& CurrentProcessInfo) override;

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
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "MonolithicDEMCoupledWallCondition" << TDim << "D";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MonolithicDEMCoupledWallCondition";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const  override{}


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
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
    }

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


    ///@}

}; // Class MonolithicDEMCoupledWallCondition


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::istream& operator >> (std::istream& rIStream,
                                  MonolithicDEMCoupledWallCondition<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MonolithicDEMCoupledWallCondition<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_MONOLITHIC_DEM_COUPLED_WALL_CONDITION_H
