//
//   Project Name:        Kratos
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_MONOLITHIC_PARTINEGRATION_PFEM2_2D_ELEM_H_INCLUDED)
#define  KRATOS_MONOLITHIC_PARTINEGRATION_PFEM2_2D_ELEM_H_INCLUDED

// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"


namespace Kratos
{

  class MonolithicAutoSlipPFEM22D
	  : public Element
   {
   public:

     /// Counted pointer of MonolithicAutoSlipPFEM22D
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(MonolithicAutoSlipPFEM22D);
    ///base type: an IndexedObject that automatically has a unique number
    ///typedef IndexedObject BaseType;
    ///Element from which it is derived
    ///typedef VMS<TDim, TNumNodes> ElementBaseType;
    ///definition of node type (default is: Node)

    //typedef Node NodeType;
    /**
     * Properties are used to store any parameters
     * related to the constitutive law
     */

    typedef Properties PropertiesType;
    ///definition of the geometry type with given NodeType
    typedef Geometry<NodeType> GeometryType;
    ///definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    //typedef typename ElementBaseType::MatrixType MatrixType;
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef std::vector<std::size_t> EquationIdVectorType;
    typedef std::vector< Dof<double>::Pointer > DofsVectorType;
    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;
	
    /// Default constructor.
     MonolithicAutoSlipPFEM22D(IndexType NewId, GeometryType::Pointer pGeometry);
     MonolithicAutoSlipPFEM22D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

     /// Destructor.
     virtual ~ MonolithicAutoSlipPFEM22D() override;


     Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

     void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

     //void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

     void AddExplicitContribution(const ProcessInfo& CurrentProcessInfo) override;

     void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

     void GetDofList(DofsVectorType& ElementalDofList,const ProcessInfo& CurrentProcessInfo) const override;

     void InitializeSolutionStep(const ProcessInfo& CurrentProcessInfo) override;



   protected:

        void CalculateLocalFractionalVelocitySystem(MatrixType& rLeftHandSideMatrix,
                                                    VectorType& rRightHandSideVector,
                                                    const ProcessInfo& rCurrentProcessInfo);

        void CalculateLocalPressureSystem(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          const ProcessInfo& rCurrentProcessInfo);

        void CalculateLocalFinalVelocitySystem(MatrixType& rLeftHandSideMatrix,
                                                VectorType& rRightHandSideVector,
                                                const ProcessInfo& rCurrentProcessInfo);

        void CalculateLocalThermalSystem(MatrixType& rLeftHandSideMatrix,
                                                VectorType& rRightHandSideVector,
                                                const ProcessInfo& rCurrentProcessInfo);

        void CalculateViscousRHS(const ProcessInfo& CurrentProcessInfo);

       	void CalculatePressureProjection(const ProcessInfo& CurrentProcessInfo);

       	void CalculateMassMatrix(const ProcessInfo& CurrentProcessInfo);

        void VelocityEquationIdVector(EquationIdVectorType& rResult,
                                      const ProcessInfo& rCurrentProcessInfo);

        void FractionalVelocityEquationIdVector(EquationIdVectorType& rResult,
                                      const ProcessInfo& rCurrentProcessInfo);

        void PressureEquationIdVector(EquationIdVectorType& rResult,
                                      const ProcessInfo& rCurrentProcessInfo);

        void ThermalEquationIdVector(EquationIdVectorType& rResult,
                                      const ProcessInfo& rCurrentProcessInfo);

        void GetVelocityDofList(DofsVectorType& rElementalDofList,
                                const ProcessInfo& rCurrentProcessInfo);

        void GetFractionalVelocityDofList(DofsVectorType& rElementalDofList,
                                const ProcessInfo& rCurrentProcessInfo);

        void GetPressureDofList(DofsVectorType& rElementalDofList,
                                const ProcessInfo& rCurrentProcessInfo);

        void GetThermalDofList(DofsVectorType& rElementalDofList,
                                const ProcessInfo& rCurrentProcessInfo);

        void AddViscousTerm(MatrixType& rDampMatrix,
                                       const BoundedMatrix<double, 3, 2>& rShapeDeriv,
                                       double Viscosity,const double Area);



       void AddViscousTerm(BoundedMatrix<double, 13, 13 > & output,
						  BoundedMatrix<double, (3), 2 >& rShapeDeriv,
						  array_1d<double,3>&  distances,
                          std::vector< Matrix >& gauss_gradients,
						  array_1d<double,3>&  viscosities,
						  array_1d<double,3>&  signs,
						  array_1d<double,3>&  volumes ,
						  const unsigned int ndivisions);

        double CalculateAirDensity(double & temperature_in_celsius);

	    double CalculateAirConductivity(double & temperature_in_celsius);


	    double CalculateAirViscosity(double & temperature_in_celsius);


	   template<class T>
	   bool InvertMatrix(const T& input, T& inverse)  ;

   private:
	friend class Serializer;

       MonolithicAutoSlipPFEM22D() : Element()
       {
       }



   }; // Class MonolithicAutoSlipPFEM22D
}  // namespace Kratos.

#endif // KRATOS_MONOLITHIC_PARTINEGRATION_PFEM2_3D_ELEM_H_INCLUDED  defined
