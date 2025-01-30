# Test line loads in staged analysis

**Authors:** [Anne van de Graaf](https://github.com/avdg81), [Gennady Markelov](https://github.com/markelov208)

## Case Specification

This test consists of a square soil domain, where X ranges from 0.0 to 5.0 m,
and Y ranges from -5.0 to 0.0 m.  It has been meshed with quadratic triangles,
which have displacement degrees of freedom as well as pore water pressure
degrees of freedom.  The left and right sides of the domain cannot move
horizontally, whereas the bottom edge has been completely fixed.  Along the top
edge, uniform vertical line loads are applied.  In the first stage, the line
load (applied to model part `First_line_load`) has a magnitude of -10.0 N/m.
In the second stage, the line load (applied to model part `Second_line_load`)
has a magnitude of -20.0 N/m.  Even though the line loads are associated with
different model parts, both of them reference the same set of nodes.  Since the
top edge is fully loaded, the total vertical reactions are expected to be equal
to 5.0 m * -10.0 N/m = -50.0 N and 5.0 m * -20.0 N/m = -100.0 N, respectively,
except for the sign, which must be reversed.

The test is performed using multi stage script and orchestrator. The ProjectParameters files are located in Legacy and Orchestrator folders for these cases, respectively.
 The multi stage script uses a ProjectParameter file for each stage when the orchestrator reads a single file that includes info for both stages (see details about the file structure below). 

````{verbatim}
{
"orchestrator" : {
"name" : "Orchestrators.KratosMultiphysics.SequentialOrchestrator",
"settings" : {
"echo_level" : 0,
"execution_list" : ["stage_1", "stage_2"],
"load_from_checkpoint" : null,
"stage_checkpoints" : true,
"stage_checkpoints_folder" : "auxiliar_files_for_python_unittest/orchestrators_files/checkpoints",
"output_validated_settings" : "test_sequential_orchestrator"
}
},
"stages" : {
"stage_1" : {
"stage_settings" : {
"analysis_stage": "Stages.KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis.GeoMechanicsAnalysis",
"problem_data": {
"problem_name":         "test",
"start_time":           0.0,
"end_time":             1.0,
"echo_level":           1,
"parallel_type":        "OpenMP",
"number_of_threads":    1
},
"solver_settings": {
<next data for the stage_1>

},
"stage_2" : {
    "stage_settings" : {
        "analysis_stage": "Stages.KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis.GeoMechanicsAnalysis",
        "problem_data": {
        "problem_name":         "test",
        "start_time":           0.0,
        "end_time":             1.0,
        "echo_level":           1,
        "parallel_type":        "OpenMP",
        "number_of_threads":    1
    },
	"problem_data": {
<next data for the stage_2>

````

It is important to move to the working folder in the test script 

       os.chdir(file_path)

before calling the orchestrator.
