/*
Kratos Multi-Physics

Copyright (c) 2015, Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
		in the documentation and/or other materials provided with the distribution.
	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
			This product includes Kratos Multi-Physics technology.
	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#if !defined(KRATOS_@{KRATOS_NAME_UPPER}_H_INCLUDED )
#define  KRATOS_@{KRATOS_NAME_UPPER}_H_INCLUDED

// System includes


// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/variables.h"

namespace Kratos
{
class @{KRATOS_NAME_CAMEL} @{KRATOS_CLASS_BASE} {
public:

  KRATOS_CLASS_POINTER_DEFINITION(@{KRATOS_NAME_CAMEL});

	typedef WeakPointerVector<Element> ParticleWeakVectorType;
	typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
	typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;

	typedef WeakPointerVector<Condition> ConditionWeakVectorType;
	typedef WeakPointerVector<Condition >::iterator ConditionWeakIteratorType;

  // Constructor void
  @{KRATOS_NAME_CAMEL}();

  // Constructor using an array of nodes
  @{KRATOS_NAME_CAMEL}( IndexType NewId, GeometryType::Pointer pGeometry );

  // Constructor using an array of nodes with properties
  @{KRATOS_NAME_CAMEL}( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

  // Destructor
  virtual ~@{KRATOS_NAME_CAMEL}();

  // Name Operations
  virtual Condition::Pointer Create(IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties ) const;


  virtual void Initialize();
  virtual void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info );
  virtual void CalculateElasticForces(VectorType& rRightHandSideVector, ProcessInfo& r_process_info );
  virtual void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& r_process_info);
  virtual void InitializeSolutionStep(ProcessInfo& r_process_infor_process_info);
  virtual void FinalizeSolutionStep(ProcessInfo& r_process_info);
  virtual void CalculateNormal(array_1d<double, 3>& rnormal);
  virtual void AddExplicitContribution(const VectorType& rRHS,
                               const Variable<VectorType>& rRHSVariable,
                               Variable<array_1d<double,3> >& rDestinationVariable,
                               const ProcessInfo& r_process_info);

  double GetYoung();
  double GetPoisson();
  double GetTgOfFrictionAngle();

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
    }


}; // class @{KRATOS_NAME_CAMEL}.

} // namespace Kratos.

#endif // KRATOS_@{KRATOS_NAME_UPPER}_H_INCLUDED  defined
