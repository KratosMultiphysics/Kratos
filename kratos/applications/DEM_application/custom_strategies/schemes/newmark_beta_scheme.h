//
// Author: Guillermo Casas gcasas@cimne.upc.edu
//

#if !defined(KRATOS_NEWMARK_BETA_SCHEME_H_INCLUDED )
#define  KRATOS_NEWMARK_BETA_SCHEME_H_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <cfloat>

// Project includes
#include "dem_integration_scheme.h"
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "DEM_application.h"

namespace Kratos {

    class NewmarkBetaScheme : public DEMIntegrationScheme {
    public:

        typedef ModelPart::NodesContainerType NodesArrayType;

        /// Pointer definition of NewmarkBetaScheme
        KRATOS_CLASS_POINTER_DEFINITION(NewmarkBetaScheme);

        /// Default constructor.
        NewmarkBetaScheme(const double gamma = 0.5, const double beta = 0.25):
            /*mGamma(gamma),*/ mBeta(beta) {}

        /// Destructor.
        virtual ~NewmarkBetaScheme() {}
        
        
        void AddSpheresVariables(ModelPart & r_model_part);
    
        void AddClustersVariables(ModelPart & r_model_part);
        
        
        void CalculateTranslationalMotion(ModelPart& model_part, NodesArrayType& pNodes, int StepFlag);

        void UpdateTranslationalVariables(
            int StepFlag,
            const Node < 3 > & i,
            array_1d<double, 3 >& coor,
            array_1d<double, 3 >& displ,
            array_1d<double, 3 >& delta_displ,
            array_1d<double, 3 >& vel,
            const array_1d<double, 3 >& initial_coor,
            const array_1d<double, 3 >& force,
            const double force_reduction_factor,
            const double mass,
            const double delta_t,
            const bool Fix_vel[3]);
        
        void UpdateRotationalVariables(
                int StepFlag,
                const Node < 3 > & i,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& angular_acceleration,
                const double delta_t,
                const bool Fix_Ang_vel[3]);        
        
        void CalculateLocalAngularAcceleration(
                                const Node < 3 > & i,
                                const double moment_of_inertia,
                                const array_1d<double, 3 >& torque, 
                                const double moment_reduction_factor,
                                array_1d<double, 3 >& angular_acceleration);
        
        void CalculateLocalAngularAccelerationByEulerEquations(
                                    const Node < 3 > & i,
                                    const array_1d<double, 3 >& local_angular_velocity,
                                    const array_1d<double, 3 >& moments_of_inertia,
                                    const array_1d<double, 3 >& local_torque, 
                                    const double moment_reduction_factor,
                                    array_1d<double, 3 >& local_angular_acceleration);
        
        /// Turn back information as a string.

        virtual std::string Info() const {
            std::stringstream buffer;
            buffer << "NewmarkBetaScheme";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const {
            rOStream << "NewmarkBetaScheme";
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const {
        }


    protected:


    private:

        //double mGamma; commented out to avoid warning
        double mBeta;


        /// Assignment operator.

        NewmarkBetaScheme& operator=(NewmarkBetaScheme const& rOther) {
            return *this;
        }

        /// Copy constructor.

        NewmarkBetaScheme(NewmarkBetaScheme const& rOther) {
            *this = rOther;
        }


        ///@}    

    }; // Class NewmarkBetaScheme


    inline std::istream& operator>>(std::istream& rIStream,
            NewmarkBetaScheme& rThis) {
        return rIStream;
    }

    inline std::ostream& operator<<(std::ostream& rOStream,
            const NewmarkBetaScheme& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

} // namespace Kratos.

#endif // KRATOS_NEWMARK_BETA_SCHEME_H_INCLUDED  defined
