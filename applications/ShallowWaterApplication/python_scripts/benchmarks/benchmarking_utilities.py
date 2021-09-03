import KratosMultiphysics as KM

from KratosMultiphysics.ShallowWaterApplication.shallow_water_analysis import ShallowWaterAnalysis

class BenchmarkingUtilities:
    def __init__(self, options = ['regular_analysis','convergence_analysis']):
        self.options = options

    @staticmethod
    def RunCase(parameters):
        model = KM.Model()
        analysis = ShallowWaterAnalysis(model, parameters)
        analysis.Run()

    @staticmethod
    def GetParametersFromListOfProcesses(processes_list, name):
        for item, list in processes_list.items():
            for process in list:
                if process['python_module'].GetString() == name:
                    return process['Parameters']

    @staticmethod
    def GetParametersFromListOfModelers(modelers_list, name):
        for modeler in modelers_list:
            if modeler['modeler_name'].GetString() == name:
                return modeler['Parameters']

    @staticmethod
    def ReplaceSettings(parameters, setting, value):
        if isinstance(value, str):
            parameters[setting].SetString(value)
        elif isinstance(value, int):
            parameters[setting].SetInt(value)
        elif isinstance(value, float):
            parameters[setting].SetDouble(value)
        else:
            raise Exception('Unsupported type')

    @classmethod
    def ReplaceProcessSettings(cls, processes_list, python_module, setting, value):
        parameters = cls.GetParametersFromListOfProcesses(processes_list, python_module)
        cls.ReplaceSettings(parameters, setting, value)

    @classmethod
    def ReplaceModelerSettings(cls, processes_list, python_module, setting, value):
        parameters = cls.GetParametersFromListOfModelers(processes_list, python_module)
        cls.ReplaceSettings(parameters, setting, value)

    def ParseArguments(self, argv, default_mode = 'regular_analysis'):
        if len(argv) > 2:
            err_msg  = 'Too many input arguments.\n'
            err_msg += self.Usage()
            raise Exception(err_msg)
        elif len(argv) == 2:
            self.mode = argv[1]
        else:
            self.mode = default_mode
            wrn_msg  = 'Setting the analysis mode to "{}"\n\n'.format(self.mode)
            wrn_msg += self.Usage()
            wrn_msg += '\n\n'
            KM.Logger.PrintWarning(wrn_msg)
    
    def Mode(self):
        return self.mode

    def Usage(self):
        usage  = 'Usage of this script:\n'
        usage += '> python MainKratos.py "mode"\n'
        usage += 'The possible modes are:\n - '
        usage += '\n - '.join(self.options)
        return usage

    def PrintUsage(self):
        KM.Logger.PrintInfo(self.Usage())

    def PrintUnknownModeMessage(self):
        msg  = 'Unknown mode. The specified value is "{}"\n'.format(self.mode)
        msg += self.Usage()
        KM.Logger.PrintInfo(msg)
