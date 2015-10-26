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
*   Last Modified by:    $Author: Feng Chun $
*   Date:                $Date: 2013-10-10 
*   Revision:            $Revision: 1.0
*
* ***********************************************************/


#if !defined(KRATOS_RIGIDFACE3D_H_INCLUDED )
#define  KRATOS_RIGIDFACE3D_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/variables.h"
#include "dem_wall.h"

namespace Kratos
{

class RigidFace3D : public DEMWall
{
public:

    // Counted pointer of RigidFace3D
    KRATOS_CLASS_POINTER_DEFINITION( RigidFace3D );
	
	
    typedef WeakPointerVector<Element> ParticleWeakVectorType; 
    typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
    typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;

    typedef WeakPointerVector<Condition> ConditionWeakVectorType; 
    typedef WeakPointerVector<Condition >::iterator ConditionWeakIteratorType;


    // Constructor void
    RigidFace3D();

    // Constructor using an array of nodes
    RigidFace3D( IndexType NewId, GeometryType::Pointer pGeometry );

    // Constructor using an array of nodes with properties
    RigidFace3D( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~RigidFace3D();


    // Name Operations

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties ) const;


    void Initialize();
    void CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo );
		
    void CalculateElasticForces(VectorType& rElasticForces, ProcessInfo& rCurrentProcessInfo);
    void CalculateNormal(array_1d<double, 3>& rnormal);
    void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo);
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);

    //double mTgOfFrictionAngle;
    //std::vector<SphericParticle*> mNeighbourSphericParticles;
    
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */



protected:
  
  
  


private:
    ///@name Static Member Variables

    /// privat variables


    // privat name Operations

	

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, DEMWall );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, DEMWall );
    }


}; // class RigidFace3D.

} // namespace Kratos.

#endif // KRATOS_RIGIDFACE3D_H_INCLUDED  defined 
