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

        double area = NodalFactor * rGeometry.DomainSize();
        // DomainSize() is the way to ask the geometry's length/area/volume (whatever is relevant for its dimension) without asking for the number of spatial dimensions first

        for(size_t itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
        {
            const NodeType& rConstNode = rGeometry[itNode];
            const double y = rConstNode.GetValue(Y_WALL); // wall distance to use in stress calculation
            if( y > 0.0  ) //&& rConstNode.GetValue(IS_STRUCTURE) != 0.0
            {
                array_1d<double,3> Vel = rGeometry[itNode].FastGetSolutionStepValue(VELOCITY);
                const array_1d<double,3>& VelMesh = rGeometry[itNode].FastGetSolutionStepValue(MESH_VELOCITY);
                Vel -= VelMesh;

                const double rho = rGeometry[itNode].FastGetSolutionStepValue(DENSITY);
                //const double nu = rGeometry[itNode].FastGetSolutionStepValue(VISCOSITY);

                double wall_vel = 0.0;
                for (size_t d = 0; d < TDim; d++)
                {
                    wall_vel += Vel[d]*Vel[d];
                }
                wall_vel = sqrt(wall_vel);

                const double slip_fac = rGeometry[itNode].FastGetSolutionStepValue(IS_SLIP);
		double edge_fac = 0.0;
		
		if(slip_fac == 30.0)
		  edge_fac = 20.0;
		else if(slip_fac == 20.0)
		  edge_fac = 1.0;
		else 
		  edge_fac = 0.01;		  
// 		double tau = rho * nu * wall_vel / y;
		double tau = edge_fac * y * rho * wall_vel * area ;		
		
//rGeometry[itNode].FastGetSolutionStepValue(SOUND_VELOCITY)=tau;
		for (size_t d = 0; d < TDim; d++)
		{
		    size_t k = itNode*BlockSize+d;
// 		    rLocalVector[k] -=  tau;
		    rLocalVector[k] -= Vel[d] * tau;
		    rLocalMatrix(k,k) += tau;
		}
                
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
