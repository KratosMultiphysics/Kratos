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

#ifndef KRATOS_LINEAR_WALL_CONDITION_H
#define KRATOS_LINEAR_WALL_CONDITION_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/process_info.h"
#include "monolithic_wall_condition.h"
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
class LinearWallCondition : public MonolithicWallCondition<TDim, TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MonolithicWallCondition
    KRATOS_CLASS_POINTER_DEFINITION(LinearWallCondition);

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
    LinearWallCondition(IndexType NewId = 0):
        MonolithicWallCondition<TDim,TNumNodes>(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    LinearWallCondition(IndexType NewId,
                        const NodesArrayType& ThisNodes):
        MonolithicWallCondition<TDim,TNumNodes>(NewId,ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    LinearWallCondition(IndexType NewId,
                        GeometryType::Pointer pGeometry):
        MonolithicWallCondition<TDim,TNumNodes>(NewId,pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    LinearWallCondition(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties):
        MonolithicWallCondition<TDim,TNumNodes>(NewId,pGeometry,pProperties)
    {
    }

    /// Copy constructor.
    LinearWallCondition(LinearWallCondition const& rOther):
        MonolithicWallCondition<TDim,TNumNodes>(rOther)
    {
    }

    /// Destructor.
    virtual ~LinearWallCondition() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    LinearWallCondition & operator=(LinearWallCondition const& rOther)
    {
        Condition::operator=(rOther);

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new LinearWallCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    virtual Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        return Condition::Pointer(new LinearWallCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
    }

    virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo)
    {
        const SizeType BlockSize = TDim + 1;
        const SizeType LocalSize = BlockSize * TNumNodes;

        if (rLeftHandSideMatrix.size1() != LocalSize)
            rLeftHandSideMatrix.resize(LocalSize,LocalSize,false);

        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize,false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);
        noalias(rRightHandSideVector) = ZeroVector(LocalSize);

        this->ApplyWallLaw(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
    }
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
        std::stringstream buffer;
        buffer << "LinearWallCondition" << TDim << "D";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "LinearWallCondition";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


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

    /// Commpute the wall stress and add corresponding terms to the system contributions.
    /**
      @param rLocalMatrix Local system matrix
      @param rLocalVector Local right hand side
      */
    void ApplyWallLaw(MatrixType& rLocalMatrix,
                      VectorType& rLocalVector,
                      ProcessInfo& rCurrentProcessInfo)
    {
        GeometryType& rGeometry = this->GetGeometry();
        const size_t BlockSize = TDim + 1;
        const double NodalFactor = 1.0 / double(TDim);

        //double& DeltaTime = rCurrentProcessInfo[DELTA_TIME];

//         double area = NodalFactor * rGeometry.DomainSize();
        // DomainSize() is the way to ask the geometry's length/area/volume (whatever is relevant for its dimension) without asking for the number of spatial dimensions first

//         for(size_t itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
//         {
//             const NodeType& rConstNode = rGeometry[itNode];
//             const double y = rConstNode.GetValue(Y_WALL); // wall distance to use in stress calculation
//             if( y > 0.0  ) //&& rConstNode.GetValue(IS_STRUCTURE) != 0.0
//             {
//                 array_1d<double,3> Vel = rGeometry[itNode].FastGetSolutionStepValue(VELOCITY);
//                 const array_1d<double,3>& VelMesh = rGeometry[itNode].FastGetSolutionStepValue(MESH_VELOCITY);
//                 Vel -= VelMesh;
//
//                 const double rho = rGeometry[itNode].FastGetSolutionStepValue(DENSITY);
//                 //const double nu = rGeometry[itNode].FastGetSolutionStepValue(VISCOSITY);
//
//                 double wall_vel = 0.0;
//                 for (size_t d = 0; d < TDim; d++)
//                 {
//                     wall_vel += Vel[d]*Vel[d];
//                 }
//                 wall_vel = sqrt(wall_vel);
//
//                 const double slip_fac = rGeometry[itNode].FastGetSolutionStepValue(IS_SLIP);
// 		double edge_fac = 0.0;
//
// 		if(slip_fac == 30.0)
// 		  edge_fac = 20.0;
// 		else if(slip_fac == 20.0)
// 		  edge_fac = 5.0;
// 		else
// 		  edge_fac = 1.0;
// // 		double tau = rho * nu * wall_vel / y;
// 		double tau = edge_fac * y  * wall_vel * area ;
//
// 		//compute bitangant direction
// 		array_1d<double,3> node_normal = rGeometry[itNode].FastGetSolutionStepValue(NORMAL);
// 		node_normal /= norm_2(node_normal);
//
// 		array_1d<double,3> bi_tangant = MathUtils<double>::CrossProduct(node_normal, Vel);
// 		array_1d<double,3> tangant;
// 		double norm_bi_tangant = norm_2(bi_tangant);
// 		if(norm_bi_tangant != 0.0){
// 		  bi_tangant =  (1.00/norm_bi_tangant) * bi_tangant;
// 		  tangant = (1.00/wall_vel) * Vel;
// 		}
// 		else {
// 		  noalias(bi_tangant) = ZeroVector(3);
// 		  noalias(tangant ) = ZeroVector(3); }
//
//
// //rGeometry[itNode].FastGetSolutionStepValue(SOUND_VELOCITY)=tau;
//
// 		size_t mm = itNode*BlockSize;
//
// 		for (size_t d = 0; d < TDim; d++)
// 		{
// 		    size_t k = itNode*BlockSize+d;
// // 		    rLocalVector[k] -=  tau;
// 		    rLocalVector[k] -= (tangant[d] + bi_tangant[d]) * tau;
// // 		    rLocalMatrix(k,k) += tau;
//
//
// 		    rLocalMatrix(k,mm) +=  tau; //(tangant[d] + bi_tangant[d]) * node_normal[0] *
// 		    rLocalMatrix(k,mm+1) +=  tau;//(tangant[d] + bi_tangant[d]) * node_normal[1] *
// 		    rLocalMatrix(k,mm+2) +=  tau;//(tangant[d] + bi_tangant[d]) * node_normal[2] *
// 		}
//
//             }
//         }
//compute normal
        if( this->GetValue(IS_STRUCTURE))
        {
            array_1d<double,3> v1, v2, An;
            v1[0] = rGeometry[1].X() - rGeometry[0].X();
            v1[1] = rGeometry[1].Y() - rGeometry[0].Y();
            v1[2] = rGeometry[1].Z() - rGeometry[0].Z();

            v2[0] = rGeometry[2].X() - rGeometry[0].X();
            v2[1] = rGeometry[2].Y() - rGeometry[0].Y();
            v2[2] = rGeometry[2].Z() - rGeometry[0].Z();



            MathUtils<double>::CrossProduct(An,v1,v2);
            An *= -0.5;
            const double area = norm_2(An);
            An /= area;
            
            //define rho as the maximum one in the element
            double rho = 0.0;
            for(size_t itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
            {
                double rho_node = rGeometry[itNode].FastGetSolutionStepValue(DENSITY);
                if(rho_node > rho) rho = rho_node;
            }

            for(size_t itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
            {
                const NodeType& rConstNode = rGeometry[itNode];
                const double y = rConstNode.GetValue(Y_WALL); // wall distance to use in stress calculation
                if( y > 0.0  ) //&& rConstNode.GetValue(IS_STRUCTURE) != 0.0
                {
                    array_1d<double,3> Vel = rConstNode.FastGetSolutionStepValue(VELOCITY);
                    const array_1d<double,3>& VelMesh = rGeometry[itNode].FastGetSolutionStepValue(MESH_VELOCITY);
                    Vel -= VelMesh;

                    const double nu = rConstNode.FastGetSolutionStepValue(VISCOSITY);

                    const double slip_fac = rConstNode.FastGetSolutionStepValue(IS_SLIP);
                    double edge_fac = 0.0;

                    if(slip_fac == 30.0)
                        edge_fac = 20.0;
                    else if(slip_fac == 20.0)
                        edge_fac = 1.0;
                    else
                        edge_fac = 0.01;
//                      double tau = area * NodalFactor  * rho * nu  / y;
                    
                     const double wall_vel = norm_2(Vel);
                     const double tau = /*edge_fac * */NodalFactor * y * rho * (1.0 + wall_vel) * area;
                    
//                     KRATOS_WATCH(tau);

                    bounded_matrix<double,3,3> block = IdentityMatrix(3,3);
                    array_1d<double,3> normal = rConstNode.FastGetSolutionStepValue(NORMAL);
                    normal /= norm_2(normal);
                    noalias(block) -= outer_prod(normal,normal);




                    const array_1d<double,3> rhs_contrib = tau * prod(block,Vel);

                    for (size_t d = 0; d < TDim; d++)
                    {
                        size_t base = itNode*BlockSize;
// 		    rLocalVector[k] -=  tau;
                        rLocalVector[base+d] -= rhs_contrib[d];
                        for(unsigned int t=0; t<TDim; t++)
                            rLocalMatrix(base+d,base+t) += tau * block(d,t);
                    }

                }
            }

//         //apply pressure boundary terms
// 	array_1d<double, TNumNodes> N;
// 	N[1] = NodalFactor;
// 	N[2] = NodalFactor;
// 	N[3] = NodalFactor;
//
// 	//compute normal
// 	array_1d<double,3> v1, v2, An;
// 	v1[0] = rGeometry[1].X() - rGeometry[0].X();
// 	v1[1] = rGeometry[1].Y() - rGeometry[0].Y();
// 	v1[2] = rGeometry[1].Z() - rGeometry[0].Z();
//
// 	v2[0] = rGeometry[2].X() - rGeometry[0].X();
// 	v2[1] = rGeometry[2].Y() - rGeometry[0].Y();
// 	v2[2] = rGeometry[2].Z() - rGeometry[0].Z();
//
// 	MathUtils<double>::CrossProduct(An,v1,v2);
//         An *= 0.5;
// 	An /= norm_2(An);
//
//         unsigned int FirstRow(0), FirstCol(0);
//         for (unsigned int i = 0; i < TNumNodes; ++i) // iterate over rows
//         {
// 	  array_1d<double,3> I_nn_n_node;
// 	  array_1d<double,3> node_normal = rGeometry[i].FastGetSolutionStepValue(NORMAL);
// 	  node_normal /= -norm_2(node_normal);
//
// 	  double dot_prod = inner_prod(node_normal,An);
// 	  I_nn_n_node = (An - dot_prod*node_normal);//dot_prod*node_normal
//
// 	  const double nodal_pr = rGeometry[i].FastGetSolutionStepValue(PRESSURE);
//
// 	    for (unsigned int m = 0; m < TDim; ++m) // iterate over v components (vx,vy[,vz])
// 	    {
//
// // 		rLocalMatrix(FirstRow + m, FirstCol + TDim) =  area * I_nn_n_node[m];
// 		rLocalMatrix(FirstRow + m, FirstCol + TDim) =  area * An[m];
// // 		rLocalVector[FirstRow + m] -= area * nodal_pr * I_nn_n_node[m];
// 		rLocalVector[FirstRow + m] -= area * nodal_pr * An[m];
//
// 	    }
// 	    FirstCol += BlockSize;
// 	    FirstRow += BlockSize;
// 	}
            //double is_cutted = 0.0;
// 	const double dist = rGeometry[0].FastGetSolutionStepValue(DISTANCE);
//         for (unsigned int i = 1; i < TNumNodes; ++i) // iterate over rows
//         {
// 	 const double other_dist =  rGeometry[i].FastGetSolutionStepValue(DISTANCE);
// 	 if ( dist*other_dist < 0.0)
// 	   is_cutted = 1.0;
// 	}



            //first point
            unsigned int FirstRow(0), FirstCol(0);
            array_1d<double,3> NN_gausse_one;
            NN_gausse_one[0] = 1.0/6.0;
            NN_gausse_one[1] = 1.0/6.0;
            NN_gausse_one[2] = 2.0/3.0;
            for (unsigned int i = 0; i < TNumNodes; ++i) // iterate over rows
            {
                double ii_flag = rGeometry[i].GetValue(IS_STRUCTURE);
                for (unsigned int j = 0; j < TNumNodes; ++j) // iterate over rows
                {
                    const double nodal_pr = rGeometry[j].FastGetSolutionStepValue(PRESSURE);
                    double jj_flag = rGeometry[j].GetValue(IS_STRUCTURE);
                    if( ii_flag != 0.0 && jj_flag != 0.0)
                    {
                        for (unsigned int m = 0; m < TDim; ++m) // iterate over v components (vx,vy[,vz])
                        {
                            rLocalMatrix(FirstRow + m, FirstCol + TDim) +=  area * An[m] * NN_gausse_one[i] * NN_gausse_one[j];
                            rLocalVector[FirstRow + m] -= area * nodal_pr * An[m] * NN_gausse_one[i] * NN_gausse_one[j];
                        }
                    }
                    FirstCol += BlockSize;

                }
                FirstCol = 0;
                FirstRow += BlockSize;
            }

            //second point
            FirstRow=0;
            FirstCol=0;
            array_1d<double,3> NN_gausse_two;
            NN_gausse_two[0] = 2.0/3.0;
            NN_gausse_two[1] = 1.0/6.0;
            NN_gausse_two[2] = 1.0/6.0;
            for (unsigned int i = 0; i < TNumNodes; ++i) // iterate over rows
            {
                double ii_flag = rGeometry[i].GetValue(IS_STRUCTURE);
                for (unsigned int j = 0; j < TNumNodes; ++j) // iterate over rows
                {
                    const double nodal_pr = rGeometry[j].FastGetSolutionStepValue(PRESSURE);
                    double jj_flag = rGeometry[j].GetValue(IS_STRUCTURE);

                    if( ii_flag != 0.0 && jj_flag != 0.0)
                    {
                        for (unsigned int m = 0; m < TDim; ++m) // iterate over v components (vx,vy[,vz])
                        {
                            rLocalMatrix(FirstRow + m, FirstCol + TDim) +=  area * An[m] * NN_gausse_two[i] * NN_gausse_two[j];
                            rLocalVector[FirstRow + m] -= area * nodal_pr * An[m] * NN_gausse_two[i] * NN_gausse_two[j];
                        }
                    }
                    FirstCol += BlockSize;

                }
                FirstCol = 0;
                FirstRow += BlockSize;
            }

            //third point
            FirstRow=0;
            FirstCol=0;
            array_1d<double,3> NN_gausse_three;
            NN_gausse_three[0] = 1.0/6.0;
            NN_gausse_three[1] = 2.0/3.0;
            NN_gausse_three[2] = 1.0/6.0;
            for (unsigned int i = 0; i < TNumNodes; ++i) // iterate over rows
            {
                double ii_flag = rGeometry[i].GetValue(IS_STRUCTURE);
                for (unsigned int j = 0; j < TNumNodes; ++j) // iterate over rows
                {
                    const double nodal_pr = rGeometry[j].FastGetSolutionStepValue(PRESSURE);
                    double jj_flag = rGeometry[j].GetValue(IS_STRUCTURE);

                    if( ii_flag != 0.0 && jj_flag != 0.0)
                    {
                        for (unsigned int m = 0; m < TDim; ++m) // iterate over v components (vx,vy[,vz])
                        {
                            rLocalMatrix(FirstRow + m, FirstCol + TDim) +=  area * An[m] * NN_gausse_three[i] * NN_gausse_three[j];
                            rLocalVector[FirstRow + m] -= area * nodal_pr * An[m] * NN_gausse_three[i] * NN_gausse_three[j];
                        }
                    }
                    FirstCol += BlockSize;

                }
                FirstCol = 0;
                FirstRow += BlockSize;
            }
        }

    }


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

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    }

    virtual void load(Serializer& rSerializer)
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

}; // Class MonolithicWallCondition


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::istream& operator >> (std::istream& rIStream,
                                  LinearWallCondition<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::ostream& operator << (std::ostream& rOStream,
                                  const LinearWallCondition<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_MONOLITHIC_WALL_CONDITION_H
