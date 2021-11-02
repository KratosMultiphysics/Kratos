import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

class ConvergencePlotter:
    '''Tool for plotting convergence graphs.
    
    The source data shall be stored in an HDF5 file and
    must follow the format of ConvergenceOutputProcess.
    '''

    def __init__(self, file_name, analysis_name='analysis_000', area=1.0):
        '''Construct the plotter, read the file and initialize the variables.

        Parameters
        ----------
        file_name : str
            The file name without extension
        analysis_name : str
            The name of the dataset inside the file
        area : float
            The area of the domain. It is used to compute the average element size
        '''
        # Read the file
        file_name = file_name + '.hdf5'
        f = h5.File(file_name)
        self.ds = f[analysis_name]

        # General data
        self.dset_time_steps = self.ds["time_step"]
        self.dset_num_elems = self.ds["num_elems"]
        self.dset_num_nodes = self.ds["num_nodes"]
        self.dset_elem_sizes = np.sqrt(area / self.dset_num_elems)

        # Initialize the filter
        self.filter = self._TrueFilter()

    def AddNewFilter(self, *, label=None, time=None):
        '''Override the current filter.'''
        self.filter = self._TrueFilter()
        self.AddFilter(label=label, time=time)

    def AddFilter(self, *, label=None, time=None):
        '''Merge the current filter with the new conditions.'''
        if label is not None:
            label_filter = [dset_label.decode() == label for dset_label in self.ds["label"]]
            self._MergeFilter(label_filter)
        if time is not None:
            time_filter = [np.abs(t - time) < dt for t, dt in zip(self.ds["time"], self.dset_time_steps)]
            self._MergeFilter(time_filter)

    def Plot(self, convergence, variable, ref_variable=None, ax=None, add_labels=True, **kwargs):
        '''Add a plot to the given axes.

        Parameters
        ----------
        convergence : str
            'spatial' or 'temporal'
        variable : str
            The name of the error variable to plot
        ref_variable : str
            If it is specified, the variable will be scaled to it
        ax : Axes
            An axes instance. Optional
        add_labels : bool
            Add labels to the plot
        kwargs
            Other arguments to pass to the plot
        '''
        # Get the data
        error, increments = self._GetData(convergence, variable, ref_variable)

        # Get the labels
        x_label, y_label = self._GetLabels(convergence, add_labels)

        # Generate the figure
        if ax is None:
            ax = plt.gca()
        if len(error) > 0:
            ax.loglog(increments, error, **kwargs)
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
        else:
            print("[WARNING]: There is no data to plot")

    def Slope(self, convergence, variable, ref_variable=None):
        '''Get the average convergence slope.

        Parameters
        ----------
        convergence : str
            'spatial' or 'temporal'
        variable : str
            The name of the error variable to compute the convergence slope
        ref_variable : str
            If it is specified, the variable will be scaled to it

        Returns
        -------
        float 
            The convergence slope
        '''
        # Get the data
        error, increments = self._GetData(convergence, variable, ref_variable)

        # Compute the slope
        if len(error) > 0:
            slopes = []
            for i in range(len(increments)-1):
                d_incr = increments[i] / increments[i+1]
                d_err = error[i] / error[i+1]
                slopes.append(d_err / d_incr)
        return sum(slopes) / len(slopes)

    def PrintLatexTable(self, variable, ref_variable=None):
        '''Print to screen the contents of a LaTeX tabular.

        Parameters
        ----------
        variable : str
            The name of the error variable to include in the tabular
        ref_variable : str
            If it is specified, the variable will be scaled to it
        '''

        # Get the data
        num_nodes = self.dset_num_nodes[self.filter]
        elem_sizes = self.dset_elem_sizes[self.filter]
        time_steps = self.dset_time_steps[self.filter]
        dset_error = self._GetDsetError(variable, ref_variable)
        error = dset_error[self.filter]

        # Print the LaTeX table
        if len(error) > 0:
            header = '$n_{nodes}$ & $\\Delta x$ & $\\Delta t$ & $e_r$ \\\\ \\hline'
            row = '{:,} & {:.2} & {:.1} & {:.2} \\\\'
            print(header)
            for num, dx, dt, err in zip(num_nodes, elem_sizes, time_steps, error):
                print(row.format(num, dx, dt, err))
        else:
            print("[WARNING]: There is no data to generate the table")

    def _GetData(self, convergence, variable, ref_variable):
        dset_error = self._GetDsetError(variable, ref_variable)
        dset_increments = self._GetDsetIncrements(convergence)
        error = dset_error[self.filter]
        increments = dset_increments[self.filter]
        return error, increments

    def _GetDsetError(self, variable, ref_variable):
        dset_error = self.ds[variable]
        if ref_variable is not None:
            reference = self.ds[ref_variable]
            dset_error = np.divide(dset_error, reference)
        return dset_error

    def _GetDsetIncrements(self, convergence):
        if convergence == "spatial":
            dset_increments = self.dset_elem_sizes
        elif convergence == "temporal":
            dset_increments = self.dset_time_steps
        else:
            raise Exception("Unknown convergence type. Invoked with '{}'. The possible values are 'spatial' or 'temporal'.")
        return dset_increments

    def _GetLabels(self, convergence, add_labels):
        x_label = None
        y_label = None
        if convergence == "spatial":
            if add_labels:
                x_label = 'mesh size, $\log(\Delta x)$'
                y_label = 'rel error norm, $\log(L_2)$'
        elif convergence == "temporal":
            if add_labels:
                x_label = 'time step, $\log(\Delta t)$'
                y_label = 'rel error norm, $\log(L_2)$'
        else:
            raise Exception("Unknown convergence type. Invoked with '{}'. The possible values are 'spatial' or 'temporal'.")
        return x_label, y_label

    def _MergeFilter(self, new_filter):
        self.filter = [a and b for a, b in zip(new_filter, self.filter)]

    def _TrueFilter(self):
        return [True] * len(self.dset_num_elems)
