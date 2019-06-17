class Logger(object):
    @staticmethod
    def PrintInfo(label, *args):
        print(label, " ".join(map(str,args)))

    @staticmethod
    def PrintWarning(label, *args):
        print(label, " ".join(map(str,args)))