class StatisticalEstimator:
    """
    Estimators of statistics of a given random variable. The estimations are
    based on a set of realisations of (and possibly some prior assumptions on)
    that random variable. The estimator is able to update its estimation when
    given new observations.
    """

    def __init__(self, **keywordArgs):
        self.data = keywordArgs.get("data", [])
        self.isDataStored = keywordArgs.get("isDataStored", False)
        self._sampleCounter = keywordArgs.get("sampleCounter", len(self.data))

    def reset(self):
        """
        Forgets about all observations and return to its initial
        state. Useful when samples are not recycled (i.e. the update
        process actually discards the previous observations).
        """
        pass

    def _addData(self, newData):
        """
        Add given value(s) to data (if self.isDataStored) and increment self._sampleCounter
        """
        if self.isDataStored:
            self.data.append(newData)
        self._sampleCounter += len(newData)
        pass

    def value(self):
        """
        Returns the current estimation.
        """

    def update(self, newData):
        """
        Update the estimation with new observations.
        """

    def sampleNumber(self):
        """
        Returns the current number of samples known to the estimator.
        """
        if self.isDataStored:
            return len(self.data)
        else:
            return self._sampleCounter
