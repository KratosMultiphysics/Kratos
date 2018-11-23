//
// Author: Guillermo Casas gcasas@cimne.upc.edu
//

#if !defined(KRATOS_HYBRID_BASHFORTH_SCHEME_H_INCLUDED )
#define  KRATOS_HYBRID_BASHFORTH_SCHEME_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <cfloat>

// Project includes
#include "../DEMApplication/custom_strategies/schemes/symplectic_euler_scheme.h"
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "../DEMApplication/custom_utilities/GeometryFunctions.h"
#include "utilities/quaternion.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) HybridBashforthScheme : public SymplecticEulerScheme {
    public:

        typedef ModelPart::NodesContainerType NodesArrayType;

        /// Pointer definition of HybridBashforthScheme
        KRATOS_CLASS_POINTER_DEFINITION(HybridBashforthScheme);

        /// Default constructor.
        HybridBashforthScheme() {
            mOldVelocity[0] = 0.0;
            mOldVelocity[1] = 0.0;
            mOldVelocity[2] = 0.0;
        }

        /// Destructor.
        virtual ~HybridBashforthScheme() {}

        DEMIntegrationScheme* CloneRaw() const override {
            DEMIntegrationScheme* cloned_scheme(new HybridBashforthScheme(*this));
            return cloned_scheme;
        }

        DEMIntegrationScheme::Pointer CloneShared() const override {
            DEMIntegrationScheme::Pointer cloned_scheme(new HybridBashforthScheme(*this));
            return cloned_scheme;
        }

    void UpdateTranslationalVariables(
        int StepFlag,
        Node < 3 > & i,
        array_1d<double, 3 >& coor,
        array_1d<double, 3 >& displ,
        array_1d<double, 3 >& delta_displ,
        array_1d<double, 3 >& vel,
        const array_1d<double, 3 >& initial_coor,
        const array_1d<double, 3 >& force,
        const double force_reduction_factor,
        const double mass,
        const double delta_t,
        const bool Fix_vel[3]) override;

        /// Turn back information as a string.

        virtual std::string Info() const override {
            std::stringstream buffer;
            buffer << "HybridBashforthScheme";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const override {
            rOStream << "SymplecticEulerScheme";
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const override {
        }


    protected:


    private:

    /// Assignment operator.

        HybridBashforthScheme& operator=(HybridBashforthScheme const& rOther) {
            return *this;
        }

        /// Copy constructor.

        HybridBashforthScheme(HybridBashforthScheme const& rOther) {
            *this = rOther;
        }

        array_1d<double, 3 > mOldVelocity;

        ///@}

    }; // Class HybridBashforthScheme

    inline std::istream& operator>>(std::istream& rIStream,
            HybridBashforthScheme& rThis) {
        return rIStream;
    }

    inline std::ostream& operator<<(std::ostream& rOStream,
            const HybridBashforthScheme& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

} // namespace Kratos.

#endif // KRATOS_HYBRID_BASHFORTH_SCHEME_H_INCLUDED  defined
