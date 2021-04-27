

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DenseLayer(settings)
    
class CallbackClass:

    def __init__(self):

        self.callback = None
        pass

    def Build(self):
        return self.callback