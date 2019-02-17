/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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

/* *********************************************************
 *
 *   Last Modified by:    $Author: rrossi $
 *   Date:                $Date: 2008-10-13 07:00:53 $
 *   Revision:            $Revision: 1.4 $
 *
 * ***********************************************************/


#if !defined(KRATOS_TOTAL_LAGRANGIAN_VELOCITY_BASED_H_INCLUDED )
#define  KRATOS_TOTAL_LAGRANGIAN_VELOCITY_BASED_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "custom_elements/total_lagrangian.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"


namespace Kratos
{

class TotalLagrangianVelocityBased
    :
public TotalLagrangian
{
public:

    // Counted pointer of TotalLagrangianVelocityBased
    KRATOS_CLASS_POINTER_DEFINITION(TotalLagrangianVelocityBased);


    // Constructor using an array of nodes
    TotalLagrangianVelocityBased(IndexType NewId, GeometryType::Pointer pGeometry);

    // Constructor using an array of nodes with properties
    TotalLagrangianVelocityBased(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    // Destructor
    virtual ~TotalLagrangianVelocityBased();


    // Name Operations

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo);

    void GetDofList(
        DofsVectorType& ElementalDofList,
        ProcessInfo& rCurrentProcessInfo);

//     void Calculate( const Variable<array_1d<double,3> >& rVariable,
//                     array_1d<double,3>& Output,
//                     const ProcessInfo& rCurrentProcessInfo);
// 
//     void Calculate( const Variable<double>& rVariable,
//                     double& Output,
//                     const ProcessInfo& rCurrentProcessInfo);

protected:

	TotalLagrangianVelocityBased() {}
private:

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;
    // A private default constructor necessary for serialization
    

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  TotalLagrangian );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer,  TotalLagrangian );
    }


}; // class KRATOS_TotalLagrangianVelocityBased_H_INCLUDED.

} // namespace Kratos.

#endif // KRATOS_TOTAL_LAGRANGIAN_VELOCITY_BASED_H_INCLUDED  defined 
