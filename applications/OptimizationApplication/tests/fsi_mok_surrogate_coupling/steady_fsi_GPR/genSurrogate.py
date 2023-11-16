import os, glob, copy, math
import numpy as np
import pandas as pd
from SU2.util import ordered_bunch
from smarty import interpolation, Log, snapshot

class genSurrogate(ordered_bunch):
    
    def __init__(self,*args,**kwarg):

        # initialize ordered bunch
        super(genSurrogate,self).__init__(*args,**kwarg)

    def snaps(self):
        self.ioDataTrain = prepData(self.fluid, self.solid, (self.ioData['train']))
        self.trainSnaps = genSnaps(copy.deepcopy(self.ioDataTrain))

        self.ioDataValid = prepData(self.fluid, self.solid, (self.ioData['valid']))
        self.validSnaps = genSnaps(copy.deepcopy(self.ioDataValid))

        self.predictSnaps = copy.deepcopy(self.validSnaps)
        
    
    def surrogates(self):
        
        self.surrogateFluid = buildSurrogate('fluid', self.trainSnaps, eval(self.interpolationMethod), self.interpolationOptionsFluid)
        self.surrogateSolid = buildSurrogate('solid', self.trainSnaps, eval(self.interpolationMethod), self.interpolationOptionsSolid)

    def predict(self):
        fluidInput = self.predictSnaps.GetAttribute('fluid')['input']
        fluidOutput = self.predictSnaps.GetAttribute('fluid')['output']
        solidInput = self.predictSnaps.GetAttribute('solid')['input']
        solidOutput = self.predictSnaps.GetAttribute('solid')['output']
        initValueDirichlet = np.array([0.0])
        for var in solidOutput:
            self.predictSnaps._Variables.update({var : initValueDirichlet})

        tol = 1e-20
        u_n = np.zeros(len(self.predictSnaps.GetValues(solidOutput)))
        u = np.zeros(len(self.predictSnaps.GetValues(solidOutput)))
        i = 0
        i_max = 100
        while 1:
            p = getPredictions((fluidInput, fluidOutput),self.predictSnaps,self.surrogateFluid)
            u_n = getPredictions((solidInput, solidOutput), self.predictSnaps,self.surrogateSolid)
            i += 1
            if i != 0:
                res = np.square(np.subtract(u_n,u)).mean()#math.sqrt(np.square(np.subtract(u_n,u)).mean())
                print('iter= {},   res= {}\n'.format(i,res))
                if res <= tol or i == i_max:
                    errorPressure = abs(self.validSnaps.GetValues(fluidOutput) - self.predictSnaps.GetValues(fluidOutput))
                    errorDisplacement = abs(self.validSnaps.GetValues(solidOutput) - self.predictSnaps.GetValues(solidOutput))
                    break
            u = u_n

    def __getattr__(self,k):
        try:
            return super(genSurrogate,self).__getattr__(k)
        except AttributeError:
            raise AttributeError('Config parameter not found')

    def __getitem__(self,k):
        try:
            return super(genSurrogate,self).__getitem__(k)
        except KeyError:
            raise KeyError('Config parameter not found: %s' % k)
#: class Config


##############################################################################
# --- MAIN FUNCTION ---
##############################################################################
#def main():
#    # Create output folder
#    if not os.path.exists(userInput["outputDir"]):
#        Log("Create output directory: {0}".format(userInput["outputDir"]))
#        os.makedirs(userInput["outputDir"])
#    else:
#        Log("Output directory: {0}".format(userInput["outputDir"]))
#
#    # Build a surrogate model.
#    Log("Build a surrogate model ...")
#    surrogate = buildSurrogate(trainingSamples,
#                               interpolationMethod=userInput["interpolationMethod"],
#                               interpolationOptions=userInput["interpolationOptions"])
#    Log("Surrogate model built.")
#
#    outputFilename = os.path.join(userInput["outputDir"], "predictions.nc")
#    Log("Predict samples and write to file {0}".format(outputFilename))
#    predictions = getPredictions(newSamples, surrogate)
#    writeOutput(predictions, outputFilename, asciiTecplotFormat=True)


###############################################################################
# --- SUPPORT FUNCTIONS ---
###############################################################################

def prepData(fluid, solid, datapath):

    dataPath = datapath

    ioData = {'fluid' : {'input' : [],
                         'output': []},
              'solid' : {'input' : [],
                         'output': []}}
    
    data = pd.read_excel(dataPath)
    data_dict = data.to_dict('list')
    nSamples = len(data_dict['DoeNo'])
    for varName in (fluid['input'] + fluid['output'] + solid['input'] + solid['output']):
        if varName not in ioData:
            varNames = [name for name in data_dict.keys() if name.__contains__(varName)]
            data = np.array(list(map(data_dict.get,varNames)))
            ioData[varName] = (varNames, data.reshape([nSamples,len(varNames)]))
            if varName in fluid['input']:
                ioData['fluid']['input'] += varNames
            if varName in fluid['output']:
                ioData['fluid']['output'] += varNames
            if varName in solid['input']:
                ioData['solid']['input'] += varNames
            if varName in solid['output']:
                ioData['solid']['output'] += varNames 
    
    return ioData


def genSnaps(ioData):

    Attributes = {}
    Attributes['fluid'] = ioData['fluid']
    Attributes['solid'] = ioData['solid']
    del ioData['fluid']
    del ioData['solid']

    Variables = {}
    for key,value in ioData.items():
        for index, varName in enumerate(value[0]):
            Variables[varName] = value[1][:,index]

    return snapshot.Snapshot(variables=Variables, attributes=Attributes)


def buildSurrogate(domain, snaps, interpolationMethod, interpolationOptions=None):
    """
    Builds a surrogate model (TPS, Gaussian, etc.)
    """
    interpolationOptions = {} if interpolationOptions is None else interpolationOptions
    options = interpolationOptions.copy()
    tuneModus = options.pop("tuneModus", None)
    tuneOptions = options.pop("tuneOptions", None)
    optimizationOptions = options.pop("optimizationOptions", None)

    input = snaps.GetValues(snaps.GetAttribute(domain)['input'], asMatrix=True)
    output = snaps.GetValues(snaps.GetAttribute(domain)['output'], asMatrix=True)

    return interpolationMethod(**options).Fit(input, output, tuneModus=tuneModus,
                                              tuneOptions=tuneOptions,
                                              optimizationOptions=optimizationOptions)


def getPredictions(ioVars,snaps, surrogate):
    """
    Predict new samples at sample sites specified by the given snapshot.
    """
    controlPoints = snaps.GetValues(ioVars[0], asMatrix=True)
    values = surrogate.GetPrediction(controlPoints)

    snaps._Variables.update(dict(zip(ioVars[1],values.T)))

    return values

def writeOutput(snaps, outputFilename, asciiTecplotFormat=False):
    """
    Writes the Snapshot to file (NetCDF format). Into ASCII-Tecplot format if
    option enabled.
    """
    snaps.WriteToFile(outputFilename)
    if asciiTecplotFormat:
        outputFilename = os.path.splitext(outputFilename)[0] + ".dat"
        variables = ["x", "y", "value"]
        snaps.WriteToTECPLOTFile(outputFilename, variables, fmt='+10.6f')

