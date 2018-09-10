//
// Author: Guillermo Casas gcasas@cimne.upc.edu
//

#if !defined(KRATOS_SYMPLECTIC_EULER_OLD_VELOCITY_SCHEME_H_INCLUDED )
#define  KRATOS_SYMPLECTIC_EULER_OLD_VELOCITY_SCHEME_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <cfloat>

// Project includes
#include "../DEM_application/custom_strategies/schemes/symplectic_euler_scheme.h"
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "../DEM_application/custom_utilities/GeometryFunctions.h"
#include "utilities/quaternion.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) SymplecticEulerOldVelocityScheme : public SymplecticEulerScheme {
    public:

        typedef ModelPart::NodesContainerType NodesArrayType;

        /// Pointer definition of SymplecticEulerOldVelocityScheme
        KRATOS_CLASS_POINTER_DEFINITION(SymplecticEulerOldVelocityScheme);

        /// Default constructor.
        SymplecticEulerOldVelocityScheme() {}

        /// Destructor.
        virtual ~SymplecticEulerOldVelocityScheme() {}

        void UpdateTranslationalVariables(
                int StepFlag,
                Node < 3 >& i,
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
            buffer << "SymplecticEulerOldVelocityScheme";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const override {
            rOStream << "SymplecticEulerOldVelocityScheme";
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const override{
        }


    protected:


    private:

    /// Assignment operator.

        SymplecticEulerOldVelocityScheme& operator=(SymplecticEulerOldVelocityScheme const& rOther) {
            return *this;
        }

        /// Copy constructor.

        SymplecticEulerOldVelocityScheme(SymplecticEulerOldVelocityScheme const& rOther) {
            *this = rOther;
        }

        ///@}

    }; // Class SymplecticEulerOldVelocityScheme

    inline std::istream& operator>>(std::istream& rIStream,
            SymplecticEulerOldVelocityScheme& rThis) {
        return rIStream;
    }

    inline std::ostream& operator<<(std::ostream& rOStream,
            const SymplecticEulerOldVelocityScheme& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

} // namespace Kratos.

#endif // KRATOS_SYMPLECTIC_EULER_OLD_VELOCITY_SCHEME_H_INCLUDED  defined
