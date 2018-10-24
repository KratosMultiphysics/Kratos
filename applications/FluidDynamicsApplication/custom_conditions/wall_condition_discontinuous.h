/*
==============================================================================
Kratos Fluid Dynamics Application
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

#ifndef KRATOS_WALL_CONDITION_DISCONTINUOUS_H
#define KRATOS_WALL_CONDITION_DISCONTINUOUS_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "wall_condition.h"
#include "includes/process_info.h"
#include "includes/deprecated_variables.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

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
class WallConditionDiscontinuous : public WallCondition<TDim,TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of WallConditionDiscontinuous
    KRATOS_CLASS_POINTER_DEFINITION(WallConditionDiscontinuous);

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

    typedef VectorMap<IndexType, DataValueContainer> SolutionStepsConditionalDataContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new condition
      */
    WallConditionDiscontinuous(IndexType NewId = 0):
        WallCondition<TDim,TNumNodes>(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    WallConditionDiscontinuous(IndexType NewId,
                               const NodesArrayType& ThisNodes):
        WallCondition<TDim,TNumNodes>(NewId,ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    WallConditionDiscontinuous(IndexType NewId,
                               GeometryType::Pointer pGeometry):
        WallCondition<TDim,TNumNodes>(NewId,pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    WallConditionDiscontinuous(IndexType NewId,
                               GeometryType::Pointer pGeometry,
                               PropertiesType::Pointer pProperties):
        WallCondition<TDim,TNumNodes>(NewId,pGeometry,pProperties)
    {
    }

    /// Copy constructor.
    WallConditionDiscontinuous(WallConditionDiscontinuous const& rOther):
        WallCondition<TDim,TNumNodes>(rOther)
    {
    }

    /// Destructor.
    ~WallConditionDiscontinuous() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Copy constructor
    WallConditionDiscontinuous & operator=(WallConditionDiscontinuous const& rOther)
    {
        Condition::operator=(rOther);

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new WallConditionDiscontinuous object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared<WallConditionDiscontinuous>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }



    Condition::Pointer Create(IndexType NewId, Condition::GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared<WallConditionDiscontinuous>(NewId, pGeom, pProperties);
    }


    /// Calculate wall stress term for all nodes with IS_STRUCTURE != 0.0
    /**
      @param rDampMatrix Left-hand side matrix
      @param rRightHandSideVector Right-hand side vector
      @param rCurrentProcessInfo ProcessInfo instance (unused)
      */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        const ProcessInfo& r_process_info = rCurrentProcessInfo;
        unsigned int step = r_process_info[FRACTIONAL_STEP];
        if ( step == 1)
        {
            // Initialize local contributions
            const SizeType LocalSize = TDim * TNumNodes;

            if (rLeftHandSideMatrix.size1() != LocalSize)
                rLeftHandSideMatrix.resize(LocalSize,LocalSize);
            if (rRightHandSideVector.size() != LocalSize)
                rRightHandSideVector.resize(LocalSize);

            noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);
            noalias(rRightHandSideVector) = ZeroVector(LocalSize);

            this->ApplyInflowCondition(rLeftHandSideMatrix,rRightHandSideVector);

            this->ApplyWallLaw(rLeftHandSideMatrix,rRightHandSideVector);
        }
        else if(step == 5)
        {
            // Initialize local contributions
            const SizeType LocalSize = TNumNodes;

            if (rLeftHandSideMatrix.size1() != LocalSize)
                rLeftHandSideMatrix.resize(LocalSize,LocalSize);
            if (rRightHandSideVector.size() != LocalSize)
                rRightHandSideVector.resize(LocalSize);

            noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);
            noalias(rRightHandSideVector) = ZeroVector(LocalSize);

            if(this->GetValue(IS_STRUCTURE) == 0.0 )
            {
                const GeometryType& rGeom = this->GetGeometry();
                const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
                const unsigned int NumGauss = IntegrationPoints.size();
                Vector GaussWeights = ZeroVector(NumGauss);

                MatrixType N = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

                array_1d<double,3> Normal;
                this->CalculateNormal(Normal); //this already contains the area
                double A = norm_2(Normal);
                Normal /= A;

                for (unsigned int g = 0; g < NumGauss; g++)
                    GaussWeights[g] = 2.0*A * IntegrationPoints[g].Weight();

                for (unsigned int g = 0; g < NumGauss; g++)
                {
                    double Weight = GaussWeights[g];

                    array_1d<double,3> vgauss = N(0,g)*rGeom[0].FastGetSolutionStepValue(VELOCITY);
                    for (unsigned int iNode = 1; iNode < TNumNodes; ++iNode)
                    {
                        vgauss += N(iNode,g)*rGeom[iNode].FastGetSolutionStepValue(VELOCITY);
                    }

                    double aux =inner_prod(Normal,vgauss);

                    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
                    {
                        rRightHandSideVector[iNode] -= N(iNode,g)*Weight*aux;
                    }
                }
            }



            // Initialize local contributions
//             const SizeType LocalSize = TNumNodes;
//
//             if (rLeftHandSideMatrix.size1() != LocalSize)
//                 rLeftHandSideMatrix.resize(LocalSize,LocalSize);
//             if (rRightHandSideVector.size() != LocalSize)
//                 rRightHandSideVector.resize(LocalSize);
//
//             noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);
//             noalias(rRightHandSideVector) = ZeroVector(LocalSize);
//
//             if(this->GetValue(IS_STRUCTURE) == 0.0 )
//             {
//                 const double N = 1.0 / static_cast<double>(TNumNodes);
//                 array_1d<double,3> rNormal;
//                 this->CalculateNormal(rNormal); //this already contains the area
//
//                 for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
//                 {
//                     const array_1d<double,3>& rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
//                     for (unsigned int d = 0; d < TDim; ++d)
//                         rRightHandSideVector[iNode] -= N*rVelocity[d]*rNormal[d];
//                 }
//             }
        }
        else
        {
            if (rLeftHandSideMatrix.size1() != 0)
                rLeftHandSideMatrix.resize(0,0,false);

            if (rRightHandSideVector.size() != 0)
                rRightHandSideVector.resize(0,false);
        }
    }


    /// Check that all data required by this condition is available and reasonable
    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

        if (Check != 0)
        {
            return Check;
        }
        else
        {
            // Check that all required variables have been registered
            if(VELOCITY.Key() == 0)
                KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check if the application was correctly registered.","");
            if(MESH_VELOCITY.Key() == 0)
                KRATOS_THROW_ERROR(std::invalid_argument,"MESH_VELOCITY Key is 0. Check if the application was correctly registered.","");
            if(NORMAL.Key() == 0)
                KRATOS_THROW_ERROR(std::invalid_argument,"NORMAL Key is 0. Check if the application was correctly registered.","")
                if(IS_STRUCTURE.Key() == 0)
                    KRATOS_THROW_ERROR(std::invalid_argument,"IS_STRUCTURE Key is 0. Check if the application was correctly registered.","");
            if(Y_WALL.Key() == 0)
                KRATOS_THROW_ERROR(std::invalid_argument,"Y_WALL Key is 0. Check if the application was correctly registered.","")

                // Checks on nodes

                // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
                for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
                {

                    if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
                        KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
                    if(this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
                        KRATOS_THROW_ERROR(std::invalid_argument,"missing MESH_VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
                    if(this->GetGeometry()[i].SolutionStepsDataHas(NORMAL) == false)
                        KRATOS_THROW_ERROR(std::invalid_argument,"missing NORMAL variable on solution step data for node ",this->GetGeometry()[i].Id());
                    if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
                            this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
                            this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
                        KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
                }

            return Check;
        }

        KRATOS_CATCH("");
    }

    /// Apply condition to prevent numerical problems due to flow into the domain in unexpected places.
    /** This condition prevents problems arising from inflow in outflow areas, typically due to vortices
     *  exiting the domain.
     * @param rLocalMatrix Local LHS matrix
     * @param rLocalVector Local RHS vector
     */
    void ApplyInflowCondition(MatrixType& rLocalMatrix,
                              VectorType& rLocalVector)
    {
        if (this->GetValue(IS_STRUCTURE) == 0.0)
        {
            const unsigned int LocalSize = TNumNodes;
            const GeometryType& rGeom = this->GetGeometry();
            const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
            const unsigned int NumGauss = IntegrationPoints.size();
            Vector GaussWeights = ZeroVector(NumGauss);

            MatrixType NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

            array_1d<double,3> Normal;
            this->CalculateNormal(Normal); //this already contains the area
            double A = std::sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]+Normal[2]*Normal[2]);
            Normal /= A;

            for (unsigned int g = 0; g < NumGauss; g++)
                GaussWeights[g] = 2.0*A * IntegrationPoints[g].Weight();

            for (unsigned int g = 0; g < NumGauss; g++)
            {
                Vector N = row(NContainer,g);
                double Weight = GaussWeights[g];

                // Neumann boundary condition
                //             for (unsigned int i = 0; i < TNumNodes; i++)
                //             {
                //                 //unsigned int row = i*LocalSize;
                //                 const NodeType& rConstNode = this->GetGeometry()[i];
                //                 if ( rConstNode.IsFixed(PRESSURE)==true)
                //                 {
                //                     const double pext = rConstNode.FastGetSolutionStepValue(PRESSURE);
                //                     for (unsigned int j = 0; j < TNumNodes; j++)
                //                     {
                //                         unsigned int row = j*LocalSize;
                //                         for (unsigned int d = 0; d < TDim;d++)
                //                             rLocalVector[row+d] -= Weight*N[j]*N[i]*pext*Normal[d];
                //                     }
                //                 }
                //             }

                // Velocity inflow correction
                array_1d<double,3> Vel = ZeroVector(3);
                double Density = 0.0;

                for (unsigned int i = 0; i < TNumNodes; i++)
                {
                    const NodeType& rConstNode = this->GetGeometry()[i];
                    Vel += N[i]*rConstNode.FastGetSolutionStepValue(VELOCITY);
                    Density += N[i]*rConstNode.FastGetSolutionStepValue(DENSITY);
                }

                double Proj = Vel[0]*Normal[0] + Vel[1]*Normal[1] + Vel[2]*Normal[2];

                if (Proj < 0)
                {
                    const double W = Weight*Density*Proj;
                    for (unsigned int i = 0; i < TNumNodes; i++)
                    {
                        double row = i*LocalSize;
                        for (unsigned int j = 0; j < TNumNodes; j++)
                        {
                            double col = j*LocalSize;
                            const array_1d<double,3>& rVel = this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY);
                            for (unsigned int d = 0; d < TDim; d++)
                            {
                                double Tij = W*N[i]*N[j];
                                rLocalMatrix(row+d,col+d) -= Tij;
                                rLocalVector[row+d] += Tij*rVel[d];
                            }
                        }
                    }
                }
            }
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
    std::string Info() const override
    {
        std::stringstream buffer;
        this->PrintInfo(buffer);
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "WallConditionDiscontinuous" << TDim << "D #" << this->Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        this->pGetGeometry()->PrintData(rOStream);
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


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        typedef WallCondition<TDim,TNumNodes> BaseType;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
    }

    void load(Serializer& rSerializer) override
    {
        typedef WallCondition<TDim,TNumNodes> BaseType;
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType );
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

}; // Class WallConditionDiscontinuous


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::istream& operator >> (std::istream& rIStream,
                                  WallConditionDiscontinuous<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::ostream& operator << (std::ostream& rOStream,
                                  const WallConditionDiscontinuous<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_WALL_CONDITION_DISCONTINUOUS_H
