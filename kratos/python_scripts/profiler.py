try:
    import line_profiler
except:
    print("Line profiler not found. Please downlad from https://pypi.org/project/line-profiler")
    exit(0)

import builtins
import KratosMultiphysics.decorator_injector as DecroatorInjector

class Profiler(object):
    def __init__(self):
        self.printed_once = False
        self.lp = line_profiler.LineProfiler()
        if 'profile' not in builtins.__dict__ or builtins.__dict__['profile'] is None:
            builtins.__dict__['profile'] = self.lp

    def ProfileAllModule(self, module):
        DecroatorInjector.InjectIntoAllModule(module, builtins.__dict__['profile'])

    def ProfileModule(self, module, expression):
        DecroatorInjector.InjectIntoModule(module, builtins.__dict__['profile'], expression)

    def PrintResults(self):
        self.printed_once = True
        self.lp.print_stats()

    def __del__(self):
        if not self.printed_once:
            self.lp.print_stats()