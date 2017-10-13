//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_NURBS_POISSON_2D_ELEM_H_INCLUDED)
#define  KRATOS_NURBS_POISSON_2D_ELEM_H_INCLUDED

// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h" 

namespace Kratos
{

  class NurbsPoisson2D
	  : public Element
   {
   public:
     
     /// Counted pointer of NurbsPoisson2D
     KRATOS_CLASS_POINTER_DEFINITION(NurbsPoisson2D);

     /// Type for shape function values container
     typedef Kratos::Vector ShapeFunctionsType;

     /// Type for a matrix containing the shape function gradients
     typedef Kratos::Matrix ShapeFunctionDerivativesType;

     /// Type for an array of shape function gradient matrices
     typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;


    /// Default constructor.
     NurbsPoisson2D(IndexType NewId, GeometryType::Pointer pGeometry);
     NurbsPoisson2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

     /// Destructor.
     virtual ~ NurbsPoisson2D();


     Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

     void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
     
     void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
     
     void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

     void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);

     void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);


      
   protected:

     virtual void CalculateGeometryData(ShapeFunctionDerivativesArrayType& rDN_DX,
                                        Matrix& rNContainer,
                                        Vector& rGaussWeights);
   
   private:
	friend class Serializer;

       NurbsPoisson2D() : Element()
       {
       }
       
       
   }; // Class NurbsPoisson2D
}  // namespace Kratos.

#endif // KRATOS_POISSON_2D_ELEM_H_INCLUDED  defined

