//verification
/*
==============================================================================
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
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: jcotela $
//   Date:                $Date: 2010-10-09 10:34:00 $
//   Revision:            $Revision: 0.1 $
//
//
#if !defined(KRATOS_MONOLITHIC_NONEWTONIAN_PFEM2_3D_ELEM_H_INCLUDED )
#define  KRATOS_MONOLITHIC_NONEWTONIAN_PFEM2_3D_ELEM_H_INCLUDED

// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "custom_elements/monolithic_2fluid_3d.h"


namespace Kratos
{

class NoNewtonianMonolithicPFEM23D : public MonolithicPFEM23D
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(NoNewtonianMonolithicPFEM23D);

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
    typedef VectorMap<IndexType, DataValueContainer> SolutionStepsElementalDataContainerType;

    /// Default constructor.

     NoNewtonianMonolithicPFEM23D(IndexType NewId = 0) : MonolithicPFEM23D(NewId)
     {
	 }

     NoNewtonianMonolithicPFEM23D(IndexType NewId, const NodesArrayType& ThisNodes) : MonolithicPFEM23D(NewId, ThisNodes)
     {
     }

     NoNewtonianMonolithicPFEM23D(IndexType NewId, GeometryType::Pointer pGeometry) : MonolithicPFEM23D(NewId, pGeometry)
     {
     }

     NoNewtonianMonolithicPFEM23D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties) : MonolithicPFEM23D(NewId, pGeometry, pProperties)
     {
     }

     /// Destructor.
     virtual ~ NoNewtonianMonolithicPFEM23D()
     {
	 }

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{
    /// Create a new element of this type
    /**
     * Returns a pointer to a new TwoFluidVMS element, created using given input
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override
    {
        return Element::Pointer(new NoNewtonianMonolithicPFEM23D(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }



protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{
    virtual void AddViscousTerm(MatrixType& rDampMatrix,
        const boost::numeric::ublas::bounded_matrix<double, 4, 3>& rShapeDeriv,
        double& Viscosity,const double Area) override;
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
    double EffectiveViscosity(double DynamicViscosity,
        double YieldStress,
        const boost::numeric::ublas::bounded_matrix<double, 3+1, 3> &rDN_DX);

    double EquivalentStrainRate(const boost::numeric::ublas::bounded_matrix<double, 3+1, 3> &rDN_DX); // TDim+1,TDim

    ///@}
    ///@name Member Variables
    ///@{
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

    NoNewtonianMonolithicPFEM23D & operator=(NoNewtonianMonolithicPFEM23D const& rOther);

    /// Copy constructor.
    NoNewtonianMonolithicPFEM23D(NoNewtonianMonolithicPFEM23D const& rOther);
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

}; // Class NoNewtonianMonolithicPFEM23D
///@}
///@name Type Definitions
///@{
///@}
///@} // Fluid Dynamics Application group
} // namespace Kratos.
#endif // NoNewtonianMonolithicPFEM23D  defined
