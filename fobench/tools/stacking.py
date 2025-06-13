import numpy as np
from obspy.signal.util import stack
from scipy.signal import hilbert


def stack_2D(data, stack_type=None):
	'''
	Co-authors: Jonas Pätzel
	Description: 
		stacks data given as 2D array with dimensions (n_signals, n_samples)
		calls obspy.signal.util.stack
		e.g. for stacking channel data
	:Params:
		- data: (type: numpy): data to stack
		- stack_type (type: String or Tuple(str, int)): type of stack
		options are: 'linear', ('pw', order) or ('root', order)
		- 'linear' stack refers to mean stack, not sum!
	:Return:
		- stacked data
	'''
	if not stack_type:
		raise ValueError('Please provide a stack type')
	return stack(data, stack_type)



def stack_3D(data, stack_type=None):
    """
    Co-authors: Jonas Pätzel
    Description:
        Stacks data given as a 3D array (n_signals, n_samples, n_windows).
        Implements linear (mean), sum, and phase-weighted stacking (PWS).
    
    :Params:
        - data: (numpy.ndarray) 3D array (n_signals, n_samples, n_windows)
        - stack_type: (str or Tuple[str, int]) Type of stack:
            - 'linear': Mean stack
            - 'sum': Summation stack
            - ('pw', order): Phase-weighted stack with given order.
    
    :Return:
        - stacked_data: (numpy.ndarray) 2D array (n_signals, n_samples).
    """
    if stack_type == 'linear':
        return np.mean(data, axis=2)

    if stack_type == 'sum': 
        return np.sum(data, axis=2)

    if isinstance(stack_type, tuple) and stack_type[0] == 'pw':
        order = stack_type[1]
        n_signals, n_samples, n_windows = data.shape

        phase_stack = np.zeros((n_signals, n_samples), dtype=np.float32)
        linear_stack = np.mean(data, axis=2)

        # Compute phase stack iteratively
        for i in range(n_windows):
            analytic_signal = hilbert(data[:, :, i], axis=1)
            phase_signal = analytic_signal / np.abs(analytic_signal)
            phase_stack += np.abs(np.mean(phase_signal, axis=0))

        phase_stack /= n_windows
        phase_stack **= order

        return linear_stack * phase_stack

    raise ValueError('Please provide a stacking method of either "linear" or ("pw", order")')