/*
==============================================================================
KratosMeshingApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


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


#ifndef KRATOS_MPI_CONNECTIVITY_PRESERVE_MODELER_H
#define KRATOS_MPI_CONNECTIVITY_PRESERVE_MODELER_H

// system includes

// kratos includes
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "modeler/modeler.h"

namespace Kratos
{

class MPIConnectivityPreserveModeler : public Modeler
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(MPIConnectivityPreserveModeler);


    /// Default Constructor
    MPIConnectivityPreserveModeler();


    /// Destructor
    ~MPIConnectivityPreserveModeler();

    virtual void GenerateModelPart(ModelPart &OriginModelPart,
                                   ModelPart &DestinationModelPart,
                                   const Element &rReferenceElement,
                                   const Condition &rReferenceBoundaryCondition) override;

private:

    /// Copy constructor
    MPIConnectivityPreserveModeler(const MPIConnectivityPreserveModeler &rOther);

};

}

#endif // KRATOS_MPI_CONNECTIVITY_PRESERVE_MODELER_H
