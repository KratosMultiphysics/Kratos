//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_MID_POINT__SCHEME_H_INCLUDED )
#define  KRATOS_MID_POINT__SCHEME_H_INCLUDED

// System includes
#include <string>
#include <iostream> 

// External includes 

// Project includes
#include "dem_integration_scheme.h"
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"


namespace Kratos {

    class MidPointScheme : public DEMIntegrationScheme {
    public:
        ///@name Type Definitions
        ///@{

        typedef ModelPart::NodesContainerType NodesArrayType;
        typedef ModelPart::ElementsContainerType ElementsArrayType;

        /// Pointer definition of MidPointScheme
        KRATOS_CLASS_POINTER_DEFINITION(MidPointScheme);

        /// Default constructor.

        MidPointScheme() {}

        /// Destructor.

        virtual ~MidPointScheme() {}
        
        DEMIntegrationScheme* CloneRaw() const override {
            DEMIntegrationScheme* cloned_scheme(new MidPointScheme(*this));
            return cloned_scheme;
        }
        
        DEMIntegrationScheme::Pointer CloneShared() const override {
            DEMIntegrationScheme::Pointer cloned_scheme(new MidPointScheme(*this));
            return cloned_scheme;
        }
        
        void SetIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose = true) const override {
            if(verbose) std::cout << "\nAssigning MidPointScheme to properties " << pProp->Id() << std::endl;
            pProp->SetValue(DEM_INTEGRATION_SCHEME_POINTER, this->CloneShared());
        }

        /*void AddSpheresVariables(ModelPart & r_model_part, bool TRotationOption)  override {

            DEMIntegrationScheme::AddSpheresVariables(r_model_part, TRotationOption);

        }

        void AddClustersVariables(ModelPart & r_model_part, bool TRotationOption)  override {

            DEMIntegrationScheme::AddClustersVariables(r_model_part, TRotationOption);

        }*/

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
                const bool Fix_vel[3])  override {

            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {
                    delta_displ[k] = delta_t * (vel [k] + (0.5 * delta_t / mass) * force[k]);
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                    vel[k] += (delta_t / mass) * force[k];
                } else {
                    delta_displ[k] = delta_t * vel[k];
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            } // dimensions                                     

        }

        void UpdateRotationalVariables(
                int StepFlag,
                const Node < 3 > & i,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& angular_acceleration,
                const double delta_t,
                const bool Fix_Ang_vel[3])  override {

            for (int k = 0; k < 3; k++) {
                if (Fix_Ang_vel[k] == false) {
                    delta_rotation[k] = delta_t * (angular_velocity[k] + (0.5 * delta_t * angular_acceleration[k]));
                    rotated_angle[k] += delta_rotation[k];
                    angular_velocity[k] += delta_t * angular_acceleration[k];
                } else {
                    delta_rotation[k] = angular_velocity[k] * delta_t;
                    rotated_angle[k] += delta_rotation[k];
                }
            }
        }

        void CalculateLocalAngularAcceleration(
                const Node < 3 > & i,
                const double moment_of_inertia,
                const array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& angular_acceleration)  override {

            for (int j = 0; j < 3; j++) {
                angular_acceleration[j] = moment_reduction_factor * torque[j] / moment_of_inertia;
            }
        }

        void CalculateLocalAngularAccelerationByEulerEquations(
                const Node < 3 > & i,
                const array_1d<double, 3 >& local_angular_velocity,
                const array_1d<double, 3 >& moments_of_inertia,
                const array_1d<double, 3 >& local_torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& local_angular_acceleration)  override {

            for (int j = 0; j < 3; j++) {
                local_angular_acceleration[j] = (local_torque[j] - (local_angular_velocity[(j + 1) % 3] * moments_of_inertia[(j + 2) % 3] * local_angular_velocity[(j + 2) % 3] - local_angular_velocity[(j + 2) % 3] * moments_of_inertia[(j + 1) % 3] * local_angular_velocity[(j + 1) % 3])) / moments_of_inertia[j];
                local_angular_acceleration[j] = local_angular_acceleration[j] * moment_reduction_factor;
            }
        }

        virtual std::string Info() const override{
            std::stringstream buffer;
            buffer << "MidPointScheme";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const override{
            rOStream << "MidPointScheme";
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const override{
        }


    protected:


    private:

        /// Assignment operator.

        MidPointScheme& operator=(MidPointScheme const& rOther) {
            return *this;
        }

        /// Copy constructor.

        MidPointScheme(MidPointScheme const& rOther) {
            *this = rOther;
        }


    };

    /// input stream function

    inline std::istream& operator>>(std::istream& rIStream,
            MidPointScheme& rThis) {
        return rIStream;
    }

    /// output stream function

    inline std::ostream& operator<<(std::ostream& rOStream,
            const MidPointScheme& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

} // namespace Kratos.

#endif // KRATOS_MID_POINT__SCHEME_H_INCLUDED  defined 
