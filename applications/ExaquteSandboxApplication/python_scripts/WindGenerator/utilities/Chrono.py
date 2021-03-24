
from time import time


class Chrono(object):

    def __init__(self, keys):
        self.keys = keys
        self.Start = {}
        self.TotalTime = {}
        for key in keys: self.TotalTime[key] = 0
    
    def start(self, key):
        self.Start[key] = time()

    def finish(self, key, show_local=False):
        LocalTime = time() - self.Start[key]
        self.TotalTime[key] += LocalTime
        if show_local:
            print("(local time) {0:9s} : {1:.3E} s".format(key, LocalTime))
        return LocalTime

    def show(self, key, normalize=None):
        if normalize is None:
            print("(total time) {0:9s}: {1:.3E} s".format(key, self.TotalTime[key]))
        else:
            print("(time/iter) {0:9s}: {1:.3E} s".format(key, self.TotalTime[key]/normalize))


    def show_all(self, normalize=None):
        for key in self.keys: self.show(key, normalize=normalize)

