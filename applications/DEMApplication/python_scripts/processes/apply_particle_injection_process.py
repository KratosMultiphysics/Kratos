import KratosMultiphysics

# Import applications
import KratosMultiphysics.DEMApplication as DEM

# Other imports

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    process_settings = settings["Parameters"]

    folder_settings = KratosMultiphysics.Parameters("""
    {
        "model_part_name"      : "please_specify_model_part_name",
        "granulometry_settings" : {
            "ParticleDiameter" : 1.0,
            "ProbabilityDistribution": 'normal',
            "StandardDeviation": 0.0,
        },
        "flow_settings" : {
            "TypeOfFlowMeasurement" : false,
            "NumberOfParticles": 100,
            "InletLimitedVelocity": 2.0,
            "InletMassFlow": 100.0,
            "DenseInletOption" : false,
        },
        "injection_settings" : {
            "InVelocityModulus": 1.0,
            "InDirectionVector": [10.0, "3*t", "x+y"],
            "VelocityDeviation": 1.0e-5,
            "injector_element_type" : 'SphericParticle',
            "injected_element_type" : 'SphericParticle',
            "Excentricity": 0.1,
            "ProbabilityDistributionOfExcentricity": 'normal',
            "StandardDeviationOfExcentricity": 0.1,
            "ClusterType": 'fromFile',
            "ClusterFilename": 'custom.clu',
        },
        "injection_interval"             : [0.0, "End"]
    }
    """)

    # Detect "End" as a tag and replace it by a large number
    if process_settings.Has("interval"):
        if process_settings["interval"][1].IsString():
            if process_settings["interval"][1].GetString() == "End":
                process_settings["interval"][1].SetDouble(1e30) # = default_settings["interval"][1]
            else:
                raise Exception("The second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())

    process_settings.AddMissingParameters(folder_settings)

    if process_settings.Has("model_part_name"):
        computing_model_part = Model[process_settings["model_part_name"].GetString()]
    else: # using default name
        computing_model_part = Model["DEM"]

    return DEM.ApplyParticleInjectionProcess(computing_model_part, process_settings)
