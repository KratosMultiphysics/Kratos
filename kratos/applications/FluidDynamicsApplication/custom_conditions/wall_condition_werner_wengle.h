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

#ifndef KRATOS_WALL_CONDITION_WERNER_WENGLE_H
#define KRATOS_WALL_CONDITION_WERNER_WENGLE_H

// System includes
#include <iostream>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/process_info.h"

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

  /// Implements a power-law wall model.
  /**
     The Werner-Wengle wall layer model in "H. Werner and H. Wengle, Large-eddy simulation
     of turbulent flow over and around a cube in a plate channel, 8th Symp. on Turbulent
     Shear Flows 19-4, 1991" is implemented. The wall shear stress is calculated
     using velocity information at the intersection point of the condition's normal
     vector with its parent element's interior face. For distributed problems the 
     MetisDivideHeterogeneousInputProcess must be used with SynchronizeConditions = true.
     This ensures the condition belongs to the same partition as its parent element.
  */
  template< unsigned int TDim, unsigned int TNumNodes = TDim >
    class WallConditionWernerWengle : public Condition
    {
    public:
    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of WallConditionWernerWengle
    KRATOS_CLASS_POINTER_DEFINITION(WallConditionWernerWengle);

    typedef Node < 3 > NodeType;

    typedef Properties PropertiesType;

    typedef Geometry<NodeType> GeometryType;

    typedef Geometry<NodeType>::PointType PointType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Geometry<NodeType>::GeometriesArrayType GeometriesArrayType;

    typedef Element::WeakPointer ElementWeakPointerType;

    typedef Element::Pointer ElementPointerType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef Kratos::Vector ShapeFunctionsType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
	@param NewId Index of the new condition
    */
    WallConditionWernerWengle(IndexType NewId = 0) : Condition(NewId), mInitializeWasPerformed(false), mpElement()
    {
    }

    /// Constructor using an array of nodes.
    /**
       @param NewId Index of the new condition
       @param ThisNodes An array containing the nodes of the new condition
    */
    WallConditionWernerWengle(IndexType NewId, const NodesArrayType& ThisNodes) 
    : Condition(NewId, ThisNodes), mInitializeWasPerformed(false), mpElement()
    {
    }

    /// Constructor using Geometry.
    /**
       @param NewId Index of the new condition
       @param pGeometry Pointer to a geometry object
    */
    WallConditionWernerWengle(IndexType NewId, GeometryType::Pointer pGeometry) 
    : Condition(NewId, pGeometry), mInitializeWasPerformed(false), mpElement()
    {
    }

    /// Constructor using Properties.
    /**
       @param NewId Index of the new condition
       @param pGeometry Pointer to a geometry object
       @param pProperties Pointer to the condition's properties
    */
    WallConditionWernerWengle(IndexType NewId, GeometryType::Pointer pGeometry,
		   PropertiesType::Pointer pProperties) 
    : Condition(NewId, pGeometry, pProperties), mInitializeWasPerformed(false), mpElement()
    {
    }

    /// Copy constructor.
    WallConditionWernerWengle(WallConditionWernerWengle const& rOther) : Condition(rOther),
    mInitializeWasPerformed(rOther.mInitializeWasPerformed), mpElement(rOther.mpElement)
    {
    }

    /// Destructor.
    virtual ~WallConditionWernerWengle(){}

    ///@}
    ///@name Operators
    ///@{

    /// Copy constructor.
    WallConditionWernerWengle& operator=(WallConditionWernerWengle const& rOther)
    {
      Condition::operator=(rOther);
      mInitializeWasPerformed = rOther.mInitializeWasPerformed;
      mpElement = rOther.mpElement;
      return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new WallConditionWernerWengle object.
    /**
       @param NewId Index of the new condition
       @param ThisNodes An array containing the nodes of the new condition
       @param pProperties Pointer to the condition's properties
    */
    virtual Condition::Pointer Create(IndexType NewId, 
				      NodesArrayType const& ThisNodes, 
				      PropertiesType::Pointer pProperties) const
    {
      return Condition::Pointer(new WallConditionWernerWengle(NewId, 
				GetGeometry().Create(ThisNodes), pProperties));
    }

    /// Find the condition's parent element.
    virtual void Initialize()
    {
      KRATOS_TRY;

      if (mInitializeWasPerformed)
	return;

      mInitializeWasPerformed = true;

      double EdgeLength;
      array_1d<double,3> Edge(3,0.0);
      GeometryType& rGeom = this->GetGeometry();
      WeakPointerVector<Element> ElementCandidates; 
      for (SizeType i = 0; i < TDim; i++) 
	{
	  WeakPointerVector<Element>& rNodeElementCandidates = rGeom[i].GetValue(NEIGHBOUR_ELEMENTS);
	  for (SizeType j = 0; j < rNodeElementCandidates.size(); j++)
	    ElementCandidates.push_back(rNodeElementCandidates(j));
	}

      std::vector<IndexType> NodeIds(TNumNodes), ElementNodeIds;

      for (SizeType i=0; i < rGeom.PointsNumber(); i++)
	NodeIds[i] = rGeom[i].Id();
      
      std::sort(NodeIds.begin(), NodeIds.end());

      for (SizeType i=0; i < ElementCandidates.size(); i++)
	{
	  GeometryType& rElementGeom = ElementCandidates[i].GetGeometry();
	  ElementNodeIds.resize(rElementGeom.PointsNumber());

	  for (SizeType j=0; j < rElementGeom.PointsNumber(); j++)
	    ElementNodeIds[j] = rElementGeom[j].Id();

	  std::sort(ElementNodeIds.begin(), ElementNodeIds.end());

	  if ( std::includes(ElementNodeIds.begin(), ElementNodeIds.end(), NodeIds.begin(), NodeIds.end()) )
	    {
	      mpElement = ElementCandidates(i);
	      
	      Edge = rElementGeom[1].Coordinates() - rElementGeom[0].Coordinates();
	      mMinEdgeLength = Edge[0]*Edge[0];
	      for (SizeType j=1; j < TDim; j++)
		mMinEdgeLength += Edge[j]*Edge[j];

	      for (SizeType j=2; j < rElementGeom.PointsNumber(); j++)
		for(SizeType k=0; k < j; k++)
		  {
		    Edge = rElementGeom[j].Coordinates() - rElementGeom[k].Coordinates();
		    EdgeLength = Edge[0]*Edge[0];
		    for (SizeType l = 1; l < TDim; l++)
		      EdgeLength += Edge[l]*Edge[l];
		    mMinEdgeLength = (EdgeLength < mMinEdgeLength) ? EdgeLength : mMinEdgeLength;
		  }
	      mMinEdgeLength = sqrt(mMinEdgeLength);
	      return;
	    }
	}

      std::cout << "error in condition -> " << this->Id() << std::endl;
      KRATOS_ERROR(std::logic_error, "Condition cannot find parent element","");
      KRATOS_CATCH("");
    }

    virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, 
				       ProcessInfo& rCurrentProcessInfo)
    {
      VectorType RHS;
      this->CalculateLocalSystem(rLeftHandSideMatrix, RHS, rCurrentProcessInfo);
    }

    /// Calculate wall stress term for all nodes with IS_STRUCTURE != 0.
    /**
       @param rDampingMatrix Left-hand side matrix
       @param rRightHandSideVector Right-hand side vector
       @param rCurrentProcessInfo ProcessInfo instance (unused)
    */
    virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
				      VectorType& rRightHandSideVector,
				      ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY;

      if (mInitializeWasPerformed == false)
	{
	  std::cout << "error in condition -> " << this->Id() << std::endl;
	  KRATOS_ERROR(std::logic_error, "Condition was not initialized","");
	}

      if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1)
	{
	  // Initialize local contributions
	  const SizeType LocalSize = TDim * TNumNodes;

	  if (rLeftHandSideMatrix.size1() != LocalSize)
	    rLeftHandSideMatrix.resize(LocalSize, LocalSize);
	  if (rRightHandSideVector.size() != LocalSize)
	    rRightHandSideVector.resize(LocalSize);

	  noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
	  noalias(rRightHandSideVector) = ZeroVector(LocalSize);
	  
	  this->ApplyWallLaw(rLeftHandSideMatrix, rRightHandSideVector);
	}
      else
	{
	  if (rLeftHandSideMatrix.size1() != 0)
	    rLeftHandSideMatrix.resize(0,0,false);
	  
	  if (rRightHandSideVector.size() != 0)
	    rRightHandSideVector.resize(0,false);
	}

      KRATOS_CATCH("");
    }

    /// Check that all data required by this condition is available and reasonable.
    virtual int Check(const ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY;
      
      int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area >= 0

      if (Check != 0)
	{
	  return Check;
	}
      else
	{
	  // Check that all required variables have been registered
	  if(VELOCITY.Key() == 0)
	    KRATOS_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check if the application was correctly registered.","");
	  if(MESH_VELOCITY.Key() == 0)
	    KRATOS_ERROR(std::invalid_argument,"MESH_VELOCITY Key is 0. Check if the application was correctly registered.","");
	  if(DENSITY.Key() == 0)
	    KRATOS_ERROR(std::invalid_argument,"DENSITY Key is 0. Check if the application was correctly registered.","");
	  if(VISCOSITY.Key() == 0)
	    KRATOS_ERROR(std::invalid_argument,"VISCOSITY Key is 0. Check if the application was correctly registered.","");
	  if(NORMAL.Key() == 0)
	    KRATOS_ERROR(std::invalid_argument,"NORMAL Key is 0. Check if the application was correctly registered.","");
	  if(IS_STRUCTURE.Key() == 0)
	    KRATOS_ERROR(std::invalid_argument,"IS_STRUCTURE Key is 0. Check if the application was correctly registered.","");

	  // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
	  for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
	    {
	      if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
		KRATOS_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
	      if(this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
		KRATOS_ERROR(std::invalid_argument,"missing MESH_VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
	      if(this->GetGeometry()[i].SolutionStepsDataHas(DENSITY) == false)
		KRATOS_ERROR(std::invalid_argument,"missing DENSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
	      if(this->GetGeometry()[i].SolutionStepsDataHas(VISCOSITY) == false)
		KRATOS_ERROR(std::invalid_argument,"missing VISCOSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
	      if(this->GetGeometry()[i].SolutionStepsDataHas(NORMAL) == false)
		KRATOS_ERROR(std::invalid_argument,"missing NORMAL variable on solution step data for node ",this->GetGeometry()[i].Id());
	      if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
		 this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
		 this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
		KRATOS_ERROR(std::invalid_argument,"missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
	    }
	  
	  return Check;
	}

      KRATOS_CATCH("");
    }

    /// Provides the global indices for each one of this condition's local rows.
    /** This determines the equation ID vector for all DOFs.
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    virtual void EquationIdVector(EquationIdVectorType& rResult,
				  ProcessInfo& rCurrentProcessInfo);

    /// Returns a list of the condition's Dofs.
    /**
     * @param ConditionDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void GetDofList(DofsVectorType& ConditionDofList,
			    ProcessInfo& CurrentProcessInfo);


    /// Returns VELOCITY_X, VELOCITY_Y, (VELOCITY_Z) for each node.
    /**
     * @param Values Vector of nodal unknowns
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    virtual void GetValuesVector(Vector& Values, int Step = 0)
    {
      const SizeType LocalSize = TDim * TNumNodes;
      unsigned int LocalIndex = 0;

      if (Values.size() != LocalSize)
	Values.resize(LocalSize, false);

      for (unsigned int i = 0; i < TNumNodes; ++i)
	{
	  array_1d<double,3>& rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
	  for (unsigned int d = 0; d < TDim; ++d)
	    Values[LocalIndex++] = rVelocity[d];
	}
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
      buffer << "WallConditionWernerWengle" << TDim << "D";
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "WallConditionWernerWengle";}
    
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

    ElementPointerType pGetElement()
    {
      return mpElement.lock();
    }

    template< class TVariableType >
    void EvaluateInPoint(TVariableType& rResult,
			 const Kratos::Variable<TVariableType>& Var,
			 const ShapeFunctionsType& rShapeFunc)
    {
      GeometryType& rGeom = this->GetGeometry();
      const SizeType NumNodes = rGeom.PointsNumber();
      
      rResult = rShapeFunc[0] * rGeom[0].FastGetSolutionStepValue(Var);
      
      for(SizeType i = 1; i < NumNodes; i++)
	{
	  rResult += rShapeFunc[i] * rGeom[i].FastGetSolutionStepValue(Var);
	}
    }

    /// Calculate input parameters to wall model.
    /** 
     * @param ywall The height of the measurement point above the wall
     * @param velt The tangential velocity vector at the measurement point
     * @param area The condition's area
     */
    void CalculateWallParameters(double& ywall, array_1d<double,3>& velt, double& area)
    {
      KRATOS_TRY;

      if (TDim != 3)
	KRATOS_ERROR(std::logic_error, "Only implemented for 3D","");

      const double tol = 1.e-8;
      double det, d, w1, w2, proj;
      array_1d<double,3> b;
      MatrixType m(3,3), minv(3,3);
      ElementPointerType pElem = pGetElement();
      const array_1d<double,3>& normal = this->GetValue(NORMAL);
      const GeometriesArrayType& faces = pElem->GetGeometry().Faces();
      const array_1d<double,3>& center = this->GetGeometry().Center();
      
      ywall = 0.0;
      area  = MathUtils<double>::Norm3(normal);
      for (SizeType i=0; i < faces.size(); i++)
	{
	  const GeometryType& rFace = faces[i];

	  // rFace[0] + w1*(rFace[1] - rFace[0]) + w2*(rFace[2] - rFace[0]) = center - d*normal
	  m(0,0) = rFace[1].X() - rFace[0].X();
	  m(1,0) = rFace[1].Y() - rFace[0].Y();
	  m(2,0) = rFace[1].Z() - rFace[0].Z();
	  m(0,1) = rFace[2].X() - rFace[0].X();
	  m(1,1) = rFace[2].Y() - rFace[0].Y();
	  m(2,1) = rFace[2].Z() - rFace[0].Z();
	  m(0,2) = normal[0];
	  m(1,2) = normal[1];
	  m(2,2) = normal[2];

	  if ( fabs(MathUtils<double>::Det3(m)) < tol * pow(mMinEdgeLength,4) )
	    continue;

	  b = center - rFace[0].Coordinates();

	  MathUtils<double>::InvertMatrix3(m, minv, det);
	  w1 = minv(0,0)*b[0] + minv(0,1)*b[1] + minv(0,2)*b[2];
	  w2 = minv(1,0)*b[0] + minv(1,1)*b[1] + minv(1,2)*b[2];
	  d  = minv(2,0)*b[0] + minv(2,1)*b[1] + minv(2,2)*b[2];
	  if (w1 >= -tol && w2 >= -tol && (w1 + w2) <= 1.+tol) // check if normal intersects this face
	    {
	      // ywall = ||d*normal|| = |d| * ||normal|| = |d| * area
	      ywall = fabs(d) * area;
	      if (ywall > tol * mMinEdgeLength) // don't count condition's face
		{
		  const array_1d<double,3> v0 = rFace[0].FastGetSolutionStepValue(VELOCITY)
		    - rFace[0].FastGetSolutionStepValue(MESH_VELOCITY);
		  const array_1d<double,3> v1 = rFace[1].FastGetSolutionStepValue(VELOCITY)
		    - rFace[1].FastGetSolutionStepValue(MESH_VELOCITY);
		  const array_1d<double,3> v2 = rFace[2].FastGetSolutionStepValue(VELOCITY)
		    - rFace[2].FastGetSolutionStepValue(MESH_VELOCITY);

		  velt[0] = w1*v1[0] + w2*v2[0] * (1. - w1 - w2)*v0[0];
		  velt[1] = w1*v1[1] + w2*v2[1] * (1. - w1 - w2)*v0[1];
		  velt[2] = w1*v1[2] + w2*v2[2] * (1. - w1 - w2)*v0[2];
		  
		  // make velocity tangent
		  proj = (velt[0]*normal[0] + velt[1]*normal[1] + velt[2]*normal[2]) / (area*area);
		  velt[0] -= proj * normal[0];
		  velt[1] -= proj * normal[1];
		  velt[2] -= proj * normal[2];

		  break;
		}
	    }
	}

      KRATOS_CATCH("");
    }

    /// Compute the wall stress and add corresponding terms to the system contributions.
    /**
       @param rLocalMatrix Local system matrix
       @param rLocalVector Local right hand side
    */
    void ApplyWallLaw(MatrixType& rLocalMatrix, VectorType& rLocalVector)
    {
      const double A = 8.3;
      const double alpha = 1./7.;
      const double tol = 1.e-8;
      const unsigned int BlockSize = TDim;
      double ywall, area, rho, nu, u, um, tmp, twall, flump;
      
      array_1d<double,3> velt;
      GeometryType& rGeometry = this->GetGeometry();
      if (   rGeometry[0].GetSolutionStepValue(Y_WALL) > 0. 
	  && rGeometry[1].GetSolutionStepValue(Y_WALL) > 0.
	  && rGeometry[2].GetSolutionStepValue(Y_WALL) > 0.)
	{
	  CalculateWallParameters(ywall, velt, area);
	  ywall = (ywall > tol * mMinEdgeLength) ? ywall : tol * mMinEdgeLength;
	  u = MathUtils<double>::Norm3(velt);
	  
	  if (u > 1.e-12)
	    {
	      const ShapeFunctionsType& N = row(this->GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_1),0);
	      EvaluateInPoint(rho, DENSITY, N);
	      EvaluateInPoint(nu, VISCOSITY, N);

	      um = nu * pow(A, 2./(1.-alpha)) / (2. * ywall);

	      // linear region
	      if (u <= um) {
		twall = 2. * rho * nu * u / ywall;
	      }
	      else { // log region
		tmp   = (1.-alpha)/2. * pow(A, (1.+alpha)/(1.-alpha)) * pow(nu/ywall, 1.+alpha);
		tmp  += (1.+alpha)/A  * pow(nu/ywall, alpha) * u;
		twall = rho * pow(tmp, 2./(1.+alpha));
	      }
	      
	      flump = (area / (double) TDim) * twall;

	      for(SizeType i=0; i < rGeometry.PointsNumber(); ++i)
		{
		  const NodeType& rNode = rGeometry[i];
		  if(rNode.GetValue(IS_STRUCTURE) != 0.0)
		    {
		      velt = rNode.FastGetSolutionStepValue(VELOCITY,1) - rNode.FastGetSolutionStepValue(MESH_VELOCITY,1);
		      tmp  = MathUtils<double>::Norm3(velt);
		      tmp  = (tmp > 1.e-12) ? tmp : 1.;
		      velt = (1. / tmp) * velt;
		      for (unsigned int d=0; d < TDim; d++)
			{
			  unsigned int k = i * BlockSize + d;
			  rLocalVector[k] -= velt[d] * flump;
			  //rLocalMatrix(k,k) += flump;
			}
		    }
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
    bool mInitializeWasPerformed;
    double mMinEdgeLength;
    ElementWeakPointerType mpElement;


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

    }; // Class WallConditionWernerWengle


  ///@}
  
  ///@name Type Definitions
  ///@{
  
  
  ///@}
  ///@name Input and output
  ///@{
  

  /// input stream function
  template< unsigned int TDim, unsigned int TNumNodes >
    inline std::istream& operator >> (std::istream& rIStream,
                                      WallConditionWernerWengle<TDim,TNumNodes>& rThis)
  {
    return rIStream;
  }

  /// output stream function
  template< unsigned int TDim, unsigned int TNumNodes >
    inline std::ostream& operator << (std::ostream& rOStream,
                                      const WallConditionWernerWengle<TDim,TNumNodes>& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    
    return rOStream;
  }
  
  ///@}
  
  ///@} addtogroup block

  
}  // namespace Kratos.

#endif // KRATOS_WALL_CONDITION_WERNER_WENGLE_H
