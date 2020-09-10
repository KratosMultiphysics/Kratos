# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3

# Import XMC methods
from xmc.tools import dynamicImport
from xmc.tools import mergeTwoListsIntoOne

class HierarchyOptimiser():
    """
    This class is in charge of deciding the best hierarchy. It is able to:
    - provide the best index set (according to its implement meaning of
            best) within adjustable constraints;
    - give the tolerance splitting associated to that choice (be it fixed or
            computed)
    - provide the sample numbers per Monte Carlo index that minimise the
            cost.

    It is assumed here that the hierarchy is optimised with respect to a single QoI.
    """

    def __init__(self,**keywordArgs):
        # Attributes
        self.indexSpace = keywordArgs.get("indexSpace",None)
        self.toleranceSplittingBounds = keywordArgs.get("toleranceSplittingBounds",None)
        self.defaultHierarchy = keywordArgs.get('defaultHierarchy',None)

        # Methods
        self.optimalIndexSet = dynamicImport(keywordArgs.get("optimalIndexSet"))
        self.optimalSampleNumbers = dynamicImport(keywordArgs.get("optimalSampleNumbers"))
        self.updateHierarchySpace = dynamicImport(keywordArgs.get("updateHierarchySpace", "xmc.tools.doNothing"))
        self.toleranceSplittingMethod = dynamicImport(keywordArgs.get("toleranceSplittingMethod", "xmc.tools.doNothing"))

        # TODO - Think of a better way to do this, improve naming
        self.varianceBlender = keywordArgs.get('varianceBlender',None)
        self.isVarianceBlended = keywordArgs.get('isVarianceBlended',False)

    # TODO This should be a static method
    def inputDictionaryTemplate(self):
        input_dict = {
            'oldHierarchy':None,
            'newSampleNumber':None,
            'errorParameters':None,
            'splittingParameter':None,
            'tolerances':None,
            'models':None,
            'parametersForModel':None,
            'estimations':None,
            'costModel':None,
            'costParameters':None,
            'costEstimations':None,
            'minimalSamplesPerIndex':5,
            'variancesForHierarchy':None,
            'defaultHierarchy':None
        }
        return input_dict

    def toleranceSplitting(self, inputDict):
        # TODO THE FOLLOWING HOLDS GOOD ONLY FOR ONE TOLERANCE AND ONE SPLITTING PARAMETER
        # DO NOT USE if more are needed
        assert(len(inputDict['tolerances'])<=1),"toleranceSplitting() cannot handle more than one tolerance for now"

        if self.toleranceSplittingBounds[0]==self.toleranceSplittingBounds[1]:
            inputDict['splittingParameter'] = self.toleranceSplittingBounds[0]
        else:
            splitting_parameter = self.toleranceSplittingMethod(inputDict)
            splitting_parameter = min(max(splitting_parameter,self.toleranceSplittingBounds[0]),self.toleranceSplittingBounds[1])
            inputDict['splittingParameter'] = splitting_parameter
        return inputDict

    def optimalHierarchy(self, inputDict):
        """
        Compute the best hierarchy by individually computing the best index
        set and best number of samples
        """
        inputDict = self.toleranceSplitting(inputDict)
        indices = self.optimalIndexSet(inputDict)
        # TODO - Think of a better way to do the following
        if self.isVarianceBlended is True:
            inputDict['variancesForHierarchy'] = self.varianceBlender.blend(indices,inputDict)
        sample_numbers = self.optimalSampleNumbers(inputDict, indices)
        sample_numbers = [max(inputDict['minimalSamplesPerIndex'],sample_numbers[i])
                for i in range(len(sample_numbers))]

        old_hierarchy = inputDict['oldHierarchy']
        new_hierarchy = []
        if (len(old_hierarchy[0]) == 0):
            new_hierarchy = [indices,sample_numbers]
        else:
            new_hierarchy = mergeTwoListsIntoOne(indices,sample_numbers)

        return new_hierarchy
