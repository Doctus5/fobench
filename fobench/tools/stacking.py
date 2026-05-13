'''
Functions for stacking arrays
'''

import numpy as np
from scipy.signal import hilbert
from obspy.signal.util import stack

def stack_2D(data: np.ndarray, stack_type: str | tuple[str, int]) -> np.ndarray:
    '''
    Stacks data given as 2D array with dimensions ``(n_signals, n_samples)``.
    Calls ``obspy.signal.util.stack``. Ideal for stacking channel data

    Parameters
    ----------
    data : np.ndarray
        Data to stack.
    stack_type : str | tuple[str, int]
        Stack type, options are:
            - ``'linear'``
            - ``('pw', order)``
            - ``('root', order)``
        linear stack refers to mean stack, not sum!

    Returns
    -------
    np.ndarray
        Stacked data.

    See also
    --------
    `obspy.signal.util.stack <https://docs.obspy.org/_modules/obspy/signal/util.html#stack>`_

    '''

    return stack(data, stack_type)

def stack_3D(data: np.ndarray, stack_type: str | tuple[str, int])-> np.ndarray:
    '''
    Stacks data given as a 3D array ``(n_signals, n_samples, n_windows)``.
    Implements linear (mean), sum, and phase-weighted stacking (PWS). Adapted
    from ``obspy.signal.util.stack``.

    Parameters
    ----------
    data : np.ndarray
        3D array (n_signals, n_samples, n_windows).
    stack_type : str | tuple[str, int]
        Type of stack, options are:
            - 'linear': mean stack
            - 'sum': summation stack
            - ('pw', order): phase-weighted stack with given order

    Raises
    ------
    ValueError
        Invalid stack_type given.

    Returns
    -------
    np.ndarray
        Stacked data.

    See also
    --------
    `obspy.signal.util.stack <https://docs.obspy.org/_modules/obspy/signal/util.html#stack>`_

    '''

    if stack_type == "linear":
        return np.mean(data, axis=2)

    elif stack_type == "sum":
        return np.sum(data, axis=2)

    elif isinstance(stack_type, tuple) and stack_type[0] == "pw":
        order = stack_type[1]
        n_signals, n_samples, n_windows = data.shape

        phase_stack = np.zeros((n_signals, n_samples), dtype=np.float32)
        linear_stack = np.mean(data, axis=2)

        # compute phase stack iteratively
        for i in range(n_windows):
            analytic_signal = hilbert(data[:, :, i], axis=1)
            phase_signal = analytic_signal / np.abs(analytic_signal)
            phase_stack += np.abs(np.mean(phase_signal, axis=0))

        phase_stack /= n_windows
        phase_stack **= order

        return linear_stack * phase_stack
    else:
        raise ValueError(
            "Please provide a stacking method of either 'linear', ('pw', order)"
            "or 'sum'")