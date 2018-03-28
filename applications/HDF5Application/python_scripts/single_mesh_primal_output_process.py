import temporal_output_process_factory


def Factory(settings, Model):
    if settings["Parameters"].Has("alpha_bossak"):
        alpha_bossak = settings["Parameters"]["alpha_bossak"].GetDouble()
        settings["Parameters"].RemoveValue("alpha_bossak")
    else:
        alpha_bossak = -0.3
    factory_helper = temporal_output_process_factory.PrimalOutputFactoryHelper(alpha_bossak)
    (temporal_output_process, model_part_output, list_of_results_output) = factory_helper.Execute(settings["Parameters"], Model)
    for results_output in list_of_results_output:
        temporal_output_process.AddOutput(results_output)
    return temporal_output_process
