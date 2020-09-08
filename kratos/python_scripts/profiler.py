import sys

try:
    import line_profiler
except:
    print("Line profiler not found. Please downlad from https://pypi.org/project/line-profiler")
    exit(0)

import builtins
import KratosMultiphysics.decorator_injector as DecroatorInjector

from contextlib import ContextDecorator

class CodeHook(object):
    def __init__(self, wrapped):
        self.__code__ = wrapped

class Profiler(object):
    def __init__(self):
        self.printed_once = False
        self.lp = line_profiler.LineProfiler()
        if 'profile' not in builtins.__dict__ or builtins.__dict__['profile'] is None:
            builtins.__dict__['profile'] = self.lp

    def __enter__(self):
        sys.setprofile(self.ProfileTracer)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        sys.setprofile(None)
        return False

    def ProfileTracer(self, frame, event, arg):
        if(event == 'call' and frame.f_code.co_name == "add_once"):
            self.lp.add_function(CodeHook(frame.f_code))
            # print(frame.f_locals)

    def ProfileAllModule(self, module):
        DecroatorInjector.InjectIntoAllContainer(module, builtins.__dict__['profile'])

    def ProfileModule(self, module, expression):
        DecroatorInjector.InjectIntoContainer(module, builtins.__dict__['profile'], expression)

    def PrintResults(self):
        self.printed_once = True
        self.lp.print_stats()

    def __del__(self):
        if not self.printed_once:
            self.lp.print_stats()

class ContextProfiler(ContextDecorator):
    def __enter__(self):
        self.lp = line_profiler.LineProfiler()
        return self.lp

    def __exit__(self, exc_type, exc_value, traceback):
        return False