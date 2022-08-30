__version__ = "2.0.0-dev"
# Alphabetical order
# TODO add method and class definition modules

from .bayesianEstimator import BayesianEstimator
from .errorEstimator import ErrorEstimator
from .estimationAssembler import EstimationAssembler
from .hierarchyOptimiser import HierarchyOptimiser
from .modelEstimator import ModelEstimator
from .momentEstimator import (
    MomentEstimator,
    CombinedMomentEstimator,
    MultiMomentEstimator,
    MultiCombinedMomentEstimator,
)
from .monoCriterion import MonoCriterion
from .monteCarloIndex import MonteCarloIndex
from .monteCarloSampler import MonteCarloSampler
from .multiCriterion import MultiCriterion
from .randomGeneratorWrapper import RandomGeneratorWrapper
from .sampleGenerator import SampleGenerator
from .classDefs_solverWrapper.singleLevelRNGSolverWrapper import SingleLevelRNGSolverWrapper
from .classDefs_solverWrapper.multiLevelRNGSolverWrapper import MultiLevelRNGSolverWrapper
from .statisticalEstimator import StatisticalEstimator
from .xmcAlgorithm import XMCAlgorithm


# Initialise ExaQUte API
from exaqute import init as exaqute_init

# Must not be called more than once
exaqute_init()
