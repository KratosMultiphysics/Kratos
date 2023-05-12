//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#if !defined(NUMERICAL_INTEGRATION_METHOD_H_INCLUDED)
#define NUMERICAL_INTEGRATION_METHOD_H_INCLUDED

// System includes

// External includes
#include "includes/define.h"
#include "includes/model_part.h"

// Project includes

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) NumericalIntegrationMethod
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(NumericalIntegrationMethod);

      // Constructor / Destructor
      NumericalIntegrationMethod();
      virtual ~NumericalIntegrationMethod();

      // Generic parameters used in integral expression functions
      struct IntegrandParams
      {
        double p1, p2, p3, p4, p5, p6, p7, p8, p9, p10;
      };

      // Public methods
      virtual void   SetNumericalIntegrationMethodInProperties (Properties::Pointer pProp, bool verbose = true) const;
      virtual double SolveIntegral                             (void);
      void           CleanParameters                           (void);

      // Public members
      double mLimMin;  // minimum limit of integration domain
      double mLimMax;  // maximum limit of integration domain
      double mCoord;   // evaluation point in integration domain
      double mTol;     // tolerance for recursive methods
      double (*mpEvalIntegrand)(NumericalIntegrationMethod*);  // pointer to function to evaluate the integrand
      IntegrandParams mParams;

      // Clone
      virtual NumericalIntegrationMethod* CloneRaw() const {
        NumericalIntegrationMethod* cloned_utl(new NumericalIntegrationMethod(*this));
        return cloned_utl;
      }

      virtual NumericalIntegrationMethod::Pointer CloneShared() const {
        NumericalIntegrationMethod::Pointer cloned_utl(new NumericalIntegrationMethod(*this));
        return cloned_utl;
      }

    private:

      // Assignment operator / Copy constructor
      NumericalIntegrationMethod& operator=(NumericalIntegrationMethod const& rOther) {return *this;}
      NumericalIntegrationMethod(NumericalIntegrationMethod const& rOther) {*this = rOther;}

      // Serializer
      friend class Serializer;

      virtual void save(Serializer& rSerializer) const {
        //rSerializer.save("MyMemberName",myMember);
      }

      virtual void load(Serializer& rSerializer) {
        //rSerializer.load("MyMemberName",myMember);
      }

  }; // Class NumericalIntegrationMethod
} // namespace Kratos

#endif // NUMERICAL_INTEGRATION_METHOD_H_INCLUDED
