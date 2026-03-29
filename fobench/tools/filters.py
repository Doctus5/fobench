'''Contains signal filtering and decimation functions'''

import warnings
from pathlib import Path

import numpy as np
import scipy.signal as signal
from scipy.signal import (cheb2ord, cheby2, convolve, get_window, iirfilter,
                          remez, medfilt2d)
from scipy.fft import fftshift, ifftshift, fftfreq, irfft2, fft2
from pyrocko.util import decimate_coeffs
try:
    from scipy.signal import sosfilt, sosfiltfilt, zpk2sos
except ImportError:
    from ._sosfilt import _sosfilt as sosfilt
    from ._sosfilt import _zpk2sos as zpk2sos
from fobench.plotting import plotting_pyqt as pyqt

def point_filter(f_type: str = None, data: np.ndarray = None, df: float = None,
                 freq: float | tuple[float, float] = None, **options) -> np.ndarray:
    '''
    Dispatcher for various filter functions

    Parameters
    ----------
    f_type : str
        The filter type.
    data : np.ndarray
        Data to filter
    df : float
        Sampling rate in Hz
    freq : float | tuple(float, float), optional
        Filter corner frequencies
    **options
        Other filter parameters passed to the individual filter functions.

    Raises
    ------
    ValueError
        Unsopported filter type.

    Returns
    -------
    result : np.ndarray
        Filtered data.

    '''

    filter_functions = {'bandpass': bandpass,
                        'bandstop': bandstop,
                        'lowpass': lowpass,
                        'highpass': highpass,
                        'median': median_filter,
                        'lp_fir': lowpass_fir,
                        'lp_cheby2': lowpass_cheby_2,
                        'afk': afk_filter,
                        'remez_fir': remez_fir}

    filter_func = filter_functions.get(f_type)
    if filter_func is None:
        raise ValueError(f'Unsupported filter type: "{f_type}", '
                         f'choose one of:\n' + '\n'.join(f"'{key}'" for key in filter_functions.keys()))

    if f_type in {'bandpass', 'bandstop', 'remez_fir'}:
        result = filter_func(data=data, df=df, freqmin=freq[0], freqmax=freq[1], **options)
    elif f_type in {'lowpass', 'highpass', 'lp_fir', 'lp_cheby2'}:
        result = filter_func(data=data, df=df, freq=freq, **options)
    elif f_type in {'median', 'afk'}:
        result = filter_func(data=data, **options)
    return result

'''
-----------------------------------------------------------------
Obspy filtering functions
modified from `here <https://docs.obspy.org/_modules/obspy/signal/filter.html>`_
# :copyright:
#     The ObsPy Development Team (devs@obspy.org)
# :license:
#     GNU Lesser General Public License, Version 3
#     (https://www.gnu.org/copyleft/lesser.html)
# :modified by:
#     Sergio Diaz (GFZ-Potsdam, sergioad@gfz-potsdam.de)
#     Jonas Pätzel (ULB Brussels, jonas.patzel[at]ulb.be)
-----------------------------------------------------------------
'''

def bandpass(data: np.ndarray, freqmin: float , freqmax: float, df: float,
             corners: int = 4, zerophase: bool = False, design: str = 'butter'):
    '''
    Butterworth-Bandpass Filter. Filters data from ``freqmin`` to ``freqmax``
    using ``corners`` corners.

    Parameters
    ----------
    data : np.ndarray
        Data to filter.
    freqmin : float
        Low corner frequency of passband.
    freqmax : float
        High corner frequency of passband.
    df : float
        Sampling rate in Hz.
    corners : int, optional
        Filter order/corners.
    zerophase : bool, optional
        If ``True`` applies the filter forward and backwards resulting in
        twice the filter order but zero phase shift in the filtered signal
    design : str, optional
        The filter design, options are:
            - Butterworth   : 'butter' (default)
            - Chebyshev I   : 'cheby1'
            - Chebyshev II  : 'cheby2'
            - Cauer/elliptic: 'ellip'
            - Bessel/Thomson: 'bessel'

    Raises
    ------
    ValueError
        One or two negative corner freqencies.
    ValueError
        ``freqmin`` larger than ``freqmax``
    ValueError
        ``freqmin`` at or above Nyquist frequency

    Returns
    -------
    np.ndarray
        Filtered signal(s).

    See also
    --------
    :func:`~scipy.signal.iirfilter` for filter design and :func:`~scipy.signal.sosfilt`/
    :func:`~scipy.signal.sosfiltfilt` for application of the filter

    '''
    if freqmin <= 0 or freqmax <= 0:
        raise ValueError(f'Corner frequencies must be positive (got {freqmin}, {freqmax}).')
    if freqmin >= freqmax:
        raise ValueError(f'freqmin ({freqmin}) must be less than freqmax ({freqmax}).')

    nyq = df/2
    low, high = freqmin/nyq, freqmax/nyq

    if low >= 1:
        raise ValueError(f'Low corner frequency ({freqmin}) is at or above Nyquist ({nyq}).')
    if high >= 1:
        warnings.warn(f'High corner frequency ({freqmax}) is at or above Nyquist ({nyq}). Applying highpass instead.')
        return highpass(data, freq=freqmin, df=df, corners=corners, zerophase=zerophase)

    z, p, k = iirfilter(corners, [low, high], btype='band', ftype=design, output='zpk')
    sos = zpk2sos(z, p, k)

    if zerophase:
        return sosfiltfilt(sos, data, axis=0)
    else:
        return sosfilt(sos, data, axis=0)

def bandstop(data: np.ndarray, freqmin: float, freqmax: float, df: float,
             corners: int = 4, zerophase: bool = False, design: str = 'butter')-> np.ndarray:
    '''
    Butterworth-Bandstop Filter. Filters data removing data between frequencies
    ``freqmin`` and ``freqmax`` using ``corners`` corners.

    Parameters
    ----------
    data : np.ndarray
        Data to filter.
    freqmin : float
        Stop band low corner frequency.
    freqmax : float
        Stop band high corner frequency.
    df : float
        Sampling rate in Hz.
    corners : int, optional
        Filter order/corners.
    zerophase : bool, optional
        If ``True`` applies the filter forward and backwards resulting in
        twice the filter order but zero phase shift in the filtered signal
    design : str, optional
        The filter design, options are:
            - Butterworth   : 'butter' (default)
            - Chebyshev I   : 'cheby1'
            - Chebyshev II  : 'cheby2'
            - Cauer/elliptic: 'ellip'
            - Bessel/Thomson: 'bessel'

    Raises
    ------
    ValueError
        One or two negative corner freqencies.
    ValueError
        ``freqmin`` at or above Nyquist frequency
    ValueError
        ``freqmin`` larger than ``freqmax``

    Returns
    -------
    np.ndarray
        Filtered signal(s).

    See also
    --------
    :func:`~scipy.signal.iirfilter` for filter design and :func:`~scipy.signal.sosfilt`/
    :func:`~scipy.signal.sosfiltfilt` for application of the filter.

    '''
    if freqmin <= 0 or freqmax <= 0:
        raise ValueError(f'Corner frequencies must be positive (got {freqmin}, {freqmax}).')
    if freqmin >= freqmax:
        raise ValueError(f'Low corner ({freqmin} Hz) must be less than '
                         f'high corner ({freqmax} Hz).')

    nyq = df/2
    low, high = freqmin/nyq, freqmax/nyq

    if low >= 1:
        raise ValueError(f'Low corner frequency ({freqmin} Hz) is at or above Nyquist ({nyq} Hz).')
    if high >= 1:
        high = 1.0
        warnings.warn(f'High corner frequency ({freqmax} Hz) is at or above '
                      f'Nyquist ({nyq} Hz). Setting Nyquist as high corner.')

    z, p, k = iirfilter(corners, [low, high], btype='bandstop', ftype=design,
                        output='zpk')
    sos = zpk2sos(z, p, k)
    if zerophase:
        return sosfiltfilt(sos, data, axis=0)
    else:
        return sosfilt(sos, data, axis=0)


def lowpass(data: np.ndarray, freq: float, df: float, corners: int = 4,
            zerophase: bool = False, design: str = 'butter')-> np.ndarray:
    '''
    Butterworth-Lowpass Filter. Filters data removing data over certain
    frequency ``freq`` using ``corners``

    Parameters
    ----------
    data : np.ndarray
        Data to filter.
    freq : float
        Filter corner frequency.
    df : float
        Sampling rate in Hz.
    corners : int, optional
        Filter order/corners
    zerophase : bool, optional
        If ``True`` applies the filter forward and backwards resulting in
        twice the filter order but zero phase shift in the filtered signal
    design : str, optional
        The filter design, options are:
            - Butterworth   : 'butter' (default)
            - Chebyshev I   : 'cheby1'
            - Chebyshev II  : 'cheby2'
            - Cauer/elliptic: 'ellip'
            - Bessel/Thomson: 'bessel'

    Raises
    ------
    ValueError
        Negative corner freqency.

    Returns
    -------
    np.ndarray
        Filtered signal(s).

    See also
    --------
    :func:`~scipy.signal.iirfilter` for filter design and :func:`~scipy.signal.sosfilt`/
    :func:`~scipy.signal.sosfiltfilt` for application of the filter.

    '''


    nyq = df/2
    f = freq / nyq
    if f <= 0:
        raise ValueError(f'Corner frequency ({freq} Hz) must be positive.')
    if f > 1:
        f = 1.0
        warnings.warn(f'Selected corner frequency ({freq} Hz) is at or above '
                      f'Nyquist ({nyq} Hz). Setting corner at Nyquist')

    z, p, k = iirfilter(corners, f, btype='lowpass', ftype='design',
                        output='zpk')
    sos = zpk2sos(z, p, k)
    if zerophase:
        return sosfiltfilt(sos, data, axis=0)
    else:
        return sosfilt(sos, data, axis=0)


def highpass(data: np.ndarray, freq: float, df: float, corners: int = 4,
             zerophase: bool = False, design: str = 'butter')-> np.ndarray:
    '''
    Butterworth-Highpass Filter. Filters data removing data below certain
    frequency ``freq`` using ``corners`` corners.

    Parameters
    ----------
    data : np.ndarray
        Data to filter.
    freq : float
        Filter corner frequency.
    df : float
        Sampling rate in Hz.
    corners : int, optional
        Filter order/corners
    zerophase : bool, optional
        If ``True`` applies the filter forward and backwards resulting in
        twice the filter order but zero phase shift in the filtered signal
    design : str, optional
        The filter design, options are:
            - Butterworth   : 'butter' (default)
            - Chebyshev I   : 'cheby1'
            - Chebyshev II  : 'cheby2'
            - Cauer/elliptic: 'ellip'
            - Bessel/Thomson: 'bessel'

    Raises
    ------
    ValueError
        Negative corner freqency.
    ValueError
        Selected corner frequency at or above Nyquist


    Returns
    -------
    np.ndarray
        Filtered signal(s).

    See also
    --------
    :func:`~scipy.signal.iirfilter` for filter design and :func:`~scipy.signal.sosfilt`/
    :func:`~scipy.signal.sosfiltfilt` for application of the filter.

    '''

    nyq = df/2
    f = freq / nyq
    if f <= 0:
        raise ValueError(f'Corner frequency ({freq} Hz) must be positive.')
    if f >= 1:
        raise ValueError('Corner frequency ({freq} Hz) is at or above '
                         f'Nyquist ({nyq} Hz).')

    z, p, k = iirfilter(corners, f, btype='highpass', ftype=design,
                        output='zpk')
    sos = zpk2sos(z, p, k)
    if zerophase:
        return sosfiltfilt(sos, data, axis=0)
    else:
        return sosfilt(sos, data, axis=0)

def remez_fir(data: np.ndarray, freqmin: float, freqmax: float, df: float)-> np.ndarray:
    '''
    Finite impulse response (FIR) filter whose transfer function minimizes
    the maximum error between the desired gain and the realized gain in the
    specified bands using the Remez exchange algorithm.

    Warning
    -------
    Obspy: "This is experimental code. Use with caution!"

    Parameters
    ----------
    data : np.ndarray
        Data to filter.
    freqmin : float
        Low corner frequency.
    freqmax : float
        High corner frequency.
    df : float
        Sampling rate in Hz.

    Returns
    -------
    np.ndarray
        Filtered data.

    See also
    --------
    For filter description see: `here <https://docs.obspy.org/_modules/obspy/signal/filter.html#remez_fir>`_

    '''

    flt = freqmin - 0.1 * freqmin # take 10% of freqmin and freqmax as "corners"
    fut = freqmax + 0.1 * freqmax
    # bandpass between freqmin and freqmax
    filt = remez(50, np.array([0, flt, freqmin, freqmax, fut, df/2-1]),
                 np.array([0, 1, 0]), Hz=df)
    return convolve(filt, data)

def lowpass_fir(data: np.ndarray, freq: float, df: float, winlen: int = 2048) -> np.ndarray:
    '''
    FIR-Lowpass Filter. Filter data by passing data only below a certain frequency.

    Warning
    -------
    Obspy: "This is experimental code. Use with caution!"

    Parameters
    ----------
    data : np.ndarray
        Data to filter.
    freq : float
        Data below this frequency pass.
    df : float
        Sampling rate in Hz.
    winlen : int, optional
        Window length for filter in samples, must be power of 2.

    Returns
    -------
    np.ndarray
        Filetered data.

    See also
    --------
    For filter description see: `here <https://docs.obspy.org/_modules/obspy/signal/filter.html#lowpass_fir>`_

    '''

    w = np.fft.fftfreq(winlen, 1/float(df))
    myfilter = np.where((abs(w) < freq), 1., 0.) # cutoff is low-pass filter
    # ideal filter
    h = np.fft.ifft(myfilter)
    beta = 11.7 # beta implies Kaiser
    myh = np.fft.fftshift(h) * get_window(beta, winlen)
    return convolve(abs(myh), data)[winlen/2:-winlen/2]


def integer_decimation(data: np.ndarray, decimation_factor: int) -> np.ndarray:
    '''
    Downsampling by simple integer decimation. The new sampling rate is initial
    sampling rate divided by ``decimation_factor``

    Warning
    -------
    Make sure that no signal is present in frequency bands above the new
    Nyquist frequency (``fs / (2 * decimation_factor)``), e.g. by applying a
    lowpass filter beforehand.

    Parameters
    ----------
    data : np.ndarray
        Data to decimate.
    decimation_factor : int
        Decimation factor.

    Raises
    ------
    TypeError
        Decimation factor not an integer.
    ValueError
        Decimation factor < 2

    Returns
    -------
    data : np.ndarray
        Decimated data of length initial length/decimation factor.

    '''

    if not isinstance(decimation_factor, int):
        raise TypeError('Decimation factor must be an integer, got {type(decimation_factor).__name__}!')
    if decimation_factor < 2:
        raise ValueError(f'Decimation factor must be 2 or greater, got {decimation_factor}.')

    return data[::decimation_factor]


def lowpass_cheby_2(data: np.ndarray, freq: float, df: float, maxorder: int = 12,
                    ba: bool = False, freq_passband: bool = False)-> np.ndarray:
    '''
    Cheby2-Lowpass Filter. Filter data by passing data only below a certain frequency.
    The main purpose of this cheby2 filter is downsampling. This method will
    iteratively design a filter, whose pass band frequency is determined dynamically,
    such that the values above the stop band frequency are lower than -96dB.

    Parameters
    ----------
    data : np.ndarray
        Data to filter.
    freq : float
        The frequency above which signals are attenuated with 95 dB.
    df : float
        Sampling rate in Hz.
    maxorder : int, optional
        Maximal order of the designed cheby2 filter
    ba : bool, optional
        If True return only the filter coefficients (b, a) instead of filtering
    freq_passband : bool, optional
        If True return additionally to the filtered data, the iteratively
        determined pass band frequency

    Returns
    -------
    tuple[float, float]
        filter coefficients (b, a)
    np.ndarray, float
        Filtered data, determined pass band frequency
    np.ndarray
        Filtered data

    '''

    nyquist = df/2
    rp, rs, order = 1, 96, 1e99 # rp - maximum ripple of passband, rs - attenuation of stopband
    ws = freq / nyquist  # stop band frequency
    wp = ws  # pass band frequency
    if ws > 1:
        ws = 1.0
        warnings.warn('Selected corner frequency is above Nyquist. ' + \
                      'Setting Nyquist as high corner.')
    while True:
        if order <= maxorder:
            break
        wp = wp*0.99
        order, wn = cheb2ord(wp, ws, rp, rs, analog=0)
    if ba:
        return cheby2(order, rs, wn, btype='low', analog=0, output='ba')
    z, p, k = cheby2(order, rs, wn, btype='low', analog=0, output='zpk')
    sos = zpk2sos(z, p, k)
    if freq_passband:
        return sosfilt(sos, data, axis=0), wp*nyquist
    return sosfilt(sos, data, axis=0)

'''
-----------------------------------------------------------------
Other filter functions
-----------------------------------------------------------------
'''

def afk_filter(data: np.ndarray, window_size: int = 32, overlap: int = 14,
    exponent: float = 0.3, normalize_power: bool = False) -> np.ndarray:
    '''
    Adaptive frequency-wavenumber filter for denoising. Adapted from pyrockos
    lightguide package

    Parameters
    ----------
    data : np.ndarray
        Data to filter.
    window_size : int, optional
        Size of window in samples
    overlap : int, optional
        Overlap in samples. The default is 14.
    exponent : float, optional
        Scaling exponent ∈ [0, 1], controls the strength of the filter.
    normalize_power : bool, optional
        Normalize amplitude spectrum.

    Raises
    ------
    ValueError
        Window size not power of two or <= 4
    ValueError
        Overlap too large
    ValueError
        Padding does not match data shape

    Returns
    -------
    np.ndarray
        Filtered data.

    See also
    --------
    `Original code <https://github.com/pyrocko/lightguide/blob/main/lightguide/afk_filter_python.py>`_
    `Isken et al. 2022 <https://doi.org/10.1093/gji/ggac229>`_

    '''

    def triangular_taper(size: int, plateau: int) -> np.ndarray:
        '''Builds triangular taper'''
        if plateau > size:
            raise ValueError('Plateau cannot be larger than size.')
        if size % 2 or plateau % 2:
            raise ValueError('Size and plateau have to be even.')
        ramp_size = int((size - plateau) / 2)
        ramp = np.linspace(0.0, 1.0, ramp_size + 2)[1:-1]
        window = np.ones(size)
        window[:ramp_size] = ramp
        window[size - ramp_size :] = ramp[::-1]
        return window * window[:, np.newaxis]

    if np.log2(window_size) % 1 or window_size < 4:
        raise ValueError(f'Window size {window_size} must be pow(2) and > 4.')
    if overlap > window_size / 2 - 1:
        raise ValueError(f'Overlap {overlap} is too large.'
                         'Maximum overlap: window_size / 2 - 1.')
    window_stride = window_size - overlap
    window_non_overlap = window_size - 2 * overlap
    npx_x, npx_y = data.shape
    nwin_x = npx_x // window_stride
    nwin_y = npx_y // window_stride
    if nwin_x % 1 or nwin_y % 1:
        raise ValueError('Padding does not match desired data shape')
    filtered_data = np.zeros_like(data)
    taper = triangular_taper(window_size, window_non_overlap)
    for iwin_x in range(nwin_x):
        px_x = iwin_x * window_stride
        slice_x = slice(px_x, px_x + window_size)
        for iwin_y in range(nwin_y):
            px_y = iwin_y * window_stride
            slice_y = slice(px_y, px_y + window_size)
            window_data = data[slice_x, slice_y]
            window_fft = np.fft.rfft2(window_data)
            power_spec = np.abs(window_fft)
            if normalize_power:
                power_spec /= power_spec.max()
            weights = power_spec**exponent
            window_fft *= weights
            window_flt = np.fft.irfft2(window_fft)
            taper_this = taper[: window_flt.shape[0], : window_flt.shape[1]]
            window_flt *= taper_this
            filtered_data[px_x : px_x + window_flt.shape[0],
                          px_y : px_y + window_flt.shape[1]] += window_flt

    return filtered_data


def lfir235()-> tuple[float, float, float]:
    '''
    :Contributors: Javier Quinteros (GFZ-Potsdam)
    Calculate filter coefficients for FIR-235 filter based on a decimation factor.

    Parameters
    ----------
    decimation_factor : int | float
        The decimation factor to be applied

    Returns
    -------
    b, a, n : float, float , float
        The filter coefficients

    '''

    filename = Path(__file__).resolve().parent.parent / 'filter_coeffs/fir235_quinteros.txt'
    b = np.loadtxt(filename)
    a = np.array([1.0])
    n = len(b) - 1
    return b, a, n


def decimate(data: np.ndarray, factor: int, f_type: str, axis: int = 0,
             n: int = None, zero_phase: bool = True):
    '''
    :Contributors: : Marius Isken (GFZ-Potsdam)
    Decimates data with appropriate LP filtering before

    Parameters
    ----------
    data : np.ndarray
        Signal(s) to decimate.
    factor : int
        Downsampling factor.
    f_type : str
        LP filter type. Options are:
            - 'iir': Chebyshev type I of order 8
            - 'fir': 30 point FIR with Hamming window
            - 'fir-remez':
            - 'fir235':
.
    axis : int, optional
        Axis along which to perform decimation.
    n : int, optional
        Filter Order.
    zero_phase : bool, optional
        Prevents phase shift in decimated signal

    Raises
    ------
    ValueError
        Unknown filter type.

    Returns
    -------
    np.ndarray
        Decimated signal(s).

    '''

    if f_type in ('iir', 'fir'):
        return signal.decimate(x=data, q=factor, n=n, ftype=f_type, axis=axis, zero_phase=zero_phase)
    elif f_type == 'fir-remez':
        b, a, _ = decimate_coeffs(q=factor, n=n, ftype='fir-remez')
    elif f_type == 'fir235':
        b, a, _ = lfir235()
    else:
        raise ValueError(f'Unknown filter type "{f_type}", options are "fir", "iir", "fir-remez" and "fir235"')

    return decimate(x=data, q=factor, ftype=signal.dlti(b, a), axis=axis, zero_phase=zero_phase)


def median_filter(data: np.ndarray, kernel_size: int | list = 3)-> np.ndarray:
    '''
    Filters data using a 2 dimensional median filter

    Parameters
    ----------
    data : np.ndarray
        Data to filter.
    kernel_size : int | list, optional
        Size of the filter kernel, must be odd. If a scalar then used as size in
        each dimension, the default corresponds to 3x3 kernel.

    Returns
    -------
    filtered data : np.ndarray.
        The filtered data

    See also
    --------
    :func:`~scipy.signal.medfilt2d`

    '''
    return medfilt2d(data, kernel_size)

def fk_filter(data: np.ndarray, dt: float, dx: float, bands: list[dict],
              propagation: str | None = None, alpha: float = 0.3,
              plot_mode:str = 'pyqt', verbose: bool = False):
    '''
    Frequency wavenumber filter

    Parameters
    ----------
    data : np.ndarray
        Data to filter.
    dt : float
        Sampling period.
    dx : float
        Channel spacing.
    bands : list[dict]
        Bands that will be retained.
    propagation : str | None, optional
        . Keep only ``'positive'`` or ``'negative'`` wavenumbers.
    alpha : float, optional
        Tukey taper parameter
    plot_mode : str
        If set to ``'pyqt'`` plot of initial wavefield, filtered wavefield, fk
        spectrum and the fk mask will be generated
    verbose : bool
        If ``True`` additionally returns the fk spectrum and fk mask

    See also
    --------
    See :func:`~fobench.tools.filters.fk_mask` for details on the bands input.

    Returns
    -------
    data_filt : np.ndarray
        Filtered data.
    data_fk : np.ndarray
        Initial data in fk domain
    mask : np.ndarray
        The fk mask

    '''
    nt, nx = data.shape

    data_fk = fftshift(fft2(data))
    f = fftshift(fftfreq(nt, d=dt))
    k = fftshift(fftfreq(nx, d=dx))

    mask = fk_mask(bands=bands, f=f, k=k, propagation=propagation, alpha=alpha)
    mask = mask.T
    data_filt = irfft2(ifftshift(data_fk * mask, axes=1)).real

    if plot_mode == 'pyqt':
        pyqt.plot_fk(wf_ini=data, wf_filt=data_filt, wf_fk=data_fk,
                    mask=mask, f=f, k=k)

    return (data_filt, data_fk, mask) if verbose else data_filt

def fk_mask(bands: list[dict], f: np.ndarray, k: np.ndarray, propagation: str = None,
            alpha: float = 0.3) -> np.ndarray:
    '''
    Builds a Tukey-tapered fk-domain mask

    Parameters
    ----------
    bands : list of dict
        The limits of the fk-domain regions that will be retained
        Each dict can contain any of:
            fmin, fmax (frequency limits)
            kmin, kmax (wavenumber limits)
            vmin, vmax (velocity limits)
    f : np.ndarray
        Frequency axis
    k : np.ndarray
        Wavenumber axis
    propagation : str | None
        Keep only ``'positive'`` or ``'negative'`` wavenumbers
    alpha : float
        Tukey taper parameter

    Raises
    ------
    ValueError
        Alpha is negative or > 1.
    ValueError
        Minimum > Maximum for any of the limits

    Returns
    -------
    mask : np.ndarray
        The fk mask of shape ``(len(k)``, ``len(f)``).

    '''

    def tukey_band(grid: np.ndarray, low: float = None, high: float = None,
                   alpha: float = 0.3)-> np.ndarray:
        '''
        Builds a tapered band-pass mask using a Tukey window for given input grid
        and limits

        Parameters
        ----------
        grid : np.ndarray
            Input grid
        low, high : float | None
            Band limits
        alpha : float
            Fraction of window used for tapering, 0 = boxcar, 1 = Hanning window.

        Returns
        -------
        mask : np.ndarray
            The fk mask.

        '''

        mask = np.zeros_like(grid, dtype=float)
        if low is None and high is None:  # no limits requested -> let everything pass
            return np.ones_like(grid, dtype=float)
        low = grid.min() if low is None else low # only one limit requested
        high = grid.max() if high is None else high


        width = high-low
        x = (grid-low) / width  # normalize to [0,1]

        if alpha <= 0: # no tapering requested
            mask[(x >= 0) & (x <= 1)] = 1
            return mask

        half = alpha / 2
        mask[(x >= half) & (x <= 1 - half)] = 1 # region not affected by taper
        in_taper = ((x >= 0) & (x < half)) | ((x > 1 - half) & (x <= 1))
        t = np.where(x < half, x / half, (1 - x) / half)
        mask[in_taper] = 0.5 * (1 - np.cos(np.pi * t[in_taper]))
        return mask

    # initialize f, k, v grids and final mask
    if alpha < 0 or alpha > 1:
        raise ValueError(f'Alpha must be between 0 and 1, got {alpha} instead!')

    f_grid = f[None, :]
    k_grid = k[:, None]
    vel_grid = np.full_like(f_grid / np.ones_like(k_grid), 10e6, dtype=float)
    np.divide(f_grid, k_grid, out=vel_grid, where=k_grid != 0)
    final_mask = np.zeros((len(k), len(f)), dtype=float)

    # loop through filter bands, create smooth masks and add to the final output
    for band in bands:
        fmin, fmax = band.get('fmin'), band.get('fmax')
        kmin, kmax = band.get('kmin'), band.get('kmax')
        vmin, vmax = band.get('vmin'), band.get('vmax')

        for param, lo, hi in [('f', fmin, fmax), ('k', kmin, kmax), ('v', vmin, vmax)]:
            if lo is not None and hi is not None and lo >= hi:
                raise ValueError(f'{param}min must be < {param}max, got {lo} >= {hi}')

        band_mask = np.ones_like(final_mask)
        band_mask = np.minimum(band_mask, tukey_band(f_grid, fmin, fmax, alpha))
        band_mask = np.minimum(band_mask, tukey_band(np.abs(k_grid), kmin, kmax, alpha))
        band_mask = np.minimum(band_mask, tukey_band(np.abs(vel_grid), vmin, vmax, alpha))

        final_mask = np.maximum(final_mask, band_mask) # merge current mask and final mask

    # add the propagation filter
    if propagation == 'positive':
        final_mask[vel_grid <= 0] = 0
    elif propagation == 'negative':
        final_mask[vel_grid >= 0] = 0

    return final_mask