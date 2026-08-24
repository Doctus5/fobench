"""Signal processing functions except filtering."""

import numpy as np
import scipy.signal as signal
import scipy.integrate as integrate
from tqdm import tqdm, trange

def hilbert(data: np.ndarray, axis: int = 0) -> np.ndarray:

    """Computes Hilbert transform for 2D array of signals.

    Parameters
    ----------
    data : np.ndarray
        2D array holding data.
    axis : int, optional
        Axis along which to operate.

    Returns
    -------
    ht : np.ndarray
        Array holding analytic signals.

    """

    return signal.hilbert(data, axis=axis)

def envelope(data: np.ndarray, axis: int = 0)-> np.ndarray:

    """Computes envelopes for 2D array of signals

    Parameters
    ----------
    data : np.ndarray
        2D array holding data.
    axis : int, optional
        Axis along which to operate.

    Returns
    -------
    env : np.ndarray
        Array holding signal envelopes.

    """

    ht = hilbert(data, axis=axis)
    env = np.sqrt(data**2 + np.real(np.conjugate(ht)*ht)) # version of Christopher Wollin
	# env = (data**2 + ht**2)**0.5 # obspy version

    return env

def integrate_signal(data: np.ndarray, dx: int, axis: int)-> np.ndarray:

    """Integrates data in space or time using the cumulative trapezoid method.

    Parameters
    ----------
    data : np.ndarray
        Data to integrate.
    dx : int
        Time or space sampling.
    axis : int
        Axis along which to apply operation.

    Returns
    -------
    np.ndarray
        Integrated and detrended data.

    """

    integ = integrate.cumulative_trapezoid(y=data, dx=dx, axis=axis, initial=0)
    return signal.detrend(integ, axis=axis)


def differentiate_signal(data: np.ndarray, method: str, axis: int, dt: int)-> np.ndarray:

    """Differentiates data in time or space.

    Parameters
    ----------
    data : np.ndarray
        Data to differentate.
    method : str
        Sets the prefered method for differentiation. can be ``"gradient"`` or ``"diff"``,
        when using ``"diff"`` data is prepended with inital value along specified axis.
    axis : int
        Axis along which to apply operation.
    dt : int
        Sampling period of data.

    Returns
    -------
    diff : np.ndarray
        Differentiated data.

    """

    if method == "gradient":
        diff = np.gradient(data, dt, axis=axis)
    elif method == "diff":
        padding = np.take(data, [0], axis=axis)
        diff = np.diff(data, axis=axis, prepend=padding) / dt

    return diff

def detrend_signal(o_signal: np.ndarray, order: int, axis:int = -1)-> np.ndarray:

    """Detrends signal(s)

    Parameters
    ----------
    o_signal : np.ndarray
        Array containing the the signal(s).
    order : int
        Order number of the fitting polynomial curve used to detrend.
    axis : int, optional
        Axis  along which to apply the detrending.

    Returns
    -------
    new_signal : np.ndarray
        Detrended signal(s).

    """

    o_signal = np.asarray(o_signal, dtype=float)

    if order == 0:
        return signal.detrend(o_signal, axis=axis, type="constant")
    if order == 1:
        return signal.detrend(o_signal, axis=axis, type="linear")

    # when order of detrned is higher than 1
    new_signal = np.empty_like(o_signal)

    if o_signal.ndim == 1:
        t = np.arange(len(o_signal))
        trend = np.polyval(np.polyfit(t, o_signal, order), t)
        new_signal = o_signal - trend
        return new_signal

    elif o_signal.ndim == 2:
        n = o_signal.shape[axis]
        t = np.arange(n)

    for i in trange(o_signal.shape[1-axis], desc="Detrending", leave=False):
        if axis == 0:
            y = o_signal[:, i]
            trend = np.polyval(np.polyfit(t, y, order), t)
            new_signal[:, i] = y - trend
        else:
            y = o_signal[i, :]
            trend = np.polyval(np.polyfit(t, y, order), t)
            new_signal[i, :] = y - trend

    return new_signal


def demean_signal(o_signal: np.ndarray, axis: int = None)-> np.ndarray:

    """Removes mean from signal(s),

    Parameters
    ----------
    o_signal : np.ndarray
        Signal(s) to demean.
    axis : int, optional
        Axis along which to apply demeaning.

    Returns
    -------
    Demeaned signal(s).

    """

    o_signal = np.asarray(o_signal)

    if o_signal.ndim == 1:
        return o_signal - np.mean(o_signal)
    elif o_signal.ndim == 2:
        return o_signal - np.mean(o_signal, axis=axis, keepdims=True)


def get_tukey_window(M: int, alpha: float, sym: bool)-> np.ndarray:

    """Returns Tukey window for filtering and tapering.

    Parameters
    ----------
    M : int
        Number of points of window.
    alpha : float
        Shape parameter, fraction of window inside tapered region, [0-1].
    sym : bool
        Returns symmetric window when ``True``.

    Returns
    -------
    np.ndarray
        Tukey window array.

    """
    return signal.windows.tukey(M, alpha, sym)

def taper_signal(data: np.ndarray, axis: int, alpha: float = 0.1, detaper: bool = False,
                 sym: bool = True) -> np.ndarray:

    """Taper signal using Tukey window.

    Parameters
    ----------
    data : np.ndarray
        Data array on which to perform tapering.
    alpha : float
        Shape parameter, fraction of window inside tapered region, [0-1].
    axis : int
        Axis along which to perform tapering.
    detaper : bool
        Option to remove taper from signal.
    sym : bool
        Symmetric window for tapering if ``True``.

    Returns
    -------
    np.ndarray
        Tapered data array.

    """

    if axis is None: #1D case
        M = data.shape[0]
        taper = get_tukey_window(M, alpha, sym=sym)
        if not detaper:
            return data * taper
        elif detaper:
            return data / taper

    M = data.shape[axis]
    taper = get_tukey_window(M, alpha, sym=True)

    taper = taper[:, None] if axis == 0 else taper[None, :]

    if not detaper:
        return np.multiply(data, taper)
    elif detaper:
        return np.divide(data, taper)

def filt_preprocess(io_signal: np.ndarray, axis:int = None, order:int = 1,
                    sym: bool = True, alpha: float = 0.1,
                    steps: tuple[bool, bool, bool] = (True, True, True)) -> np.ndarray:

    """Performs standard pre-processing of the signal before filtering.
    The steps are: detrend, demean, and taper with Tukey window.

    Parameters
    ----------
    io_signal : np.ndarray
        Original input signal.
    axis : int, optional
        Axis along which to apply operation(s).
    order : int, optional
        Order number of the fitting polynomial curve used to detrend.
    sym : bool, optional
        Symmetric window for tapering if ``True``.
    alpha : float, optional
        Shape parameter, fraction of window inside tapered region, [0-1].
    steps : tuple[bool, bool, bool], optional
        Defines which steps to perform.

    Returns
    -------
    io_signal : np.ndarray
        Pre-processed signals.

    """

    do_detrend, do_demean, do_taper = steps

    if do_detrend:
        io_signal = detrend_signal(io_signal, order, axis=axis)
    if do_demean:
        io_signal = demean_signal(io_signal, axis=axis)
    if do_taper:
        io_signal = taper_signal(io_signal, alpha=alpha, axis=axis)

    return io_signal

def normalize_signal(data: np.ndarray, method:str = "absolute max",
                     axis: int = None, ram_window: int = None, fs: int = None,
                     total_channels: int = None, num_points: int = None)-> np.ndarray:
    """Normalizes signal(s)

    Parameters
    ----------
    data : np.ndarray
        Data array to normalize.
    method : str, optional
        Normalization method, options are:
			- ``"absolute max"``: with respect to the whole record (default)
			- ``"trace max"``: for each channel/timestep individually
			- ``"running mean"``: running absolute mean normalization
			- ``"1bit"``: 1-bit normalization.
    axis : int, optional
        Axis along which to perform operation.
    ram_window : int, optional
        Window length in seconds, only for running absolute mean normalization.
    fs : int, optional
        Sampling frequency of signal.
    total_channels : int, optional
        Total number of channels, only for running absolute mean normalization.
    num_points : int, optional
        Total number of time samples, only for running absolute mean normalization.

    Raises
    ------
    TypeError
        Running mean normalization chosen without window_length given.
    ValueError
        No valid normalization method chosen.

    Returns
    -------
    normalized_data : np.ndarray
        Normalized data array.

    """

    if method == "absolute max":
        normalized_data = (data - data.min()) / (data.max() - data.min())
        normalized_data = normalized_data * 2 -1

    elif method == "trace max":
        channel_min = data.min(axis=axis, keepdims=True)
        channel_max = data.max(axis=axis, keepdims=True)
        normalized_data = (data - channel_min)/(channel_max - channel_min)
        normalized_data = normalized_data * 2 -1

    elif method == "running mean":
        if ram_window is None:
            raise TypeError("⚠️ Please provide a window length for the running absolute mean normalization")

        normalized_data = data.copy()
        # from scipy.ndimage import uniform_filter1d

        w_len = int(fs * ram_window)
        t_axis = 0 if axis == 1 else 1

        # weight = uniform_filter1d(np.abs(data), size=w_len, axis=t_axis, mode="nearest")
        # weight[weight == 0] = 1
        # normalized_data = data / weight

        for i in tqdm(range(total_channels), desc="Running mean normalization", leave=False):
            for segment_start in range(num_points - w_len + 1):
                segment_end = segment_start + w_len

                if t_axis == 0:
                    segment = normalized_data[segment_start:segment_end, i]
                    weight = np.mean(np.abs(segment))
                    normalized_data[segment_start:segment_end, i] /= weight
                else:
                    segment = normalized_data[i,segment_start:segment_end]
                    weight = np.mean(np.abs(segment))
                    normalized_data[i, segment_start:segment_end] /= weight

    elif method == "1bit":
        normalized_data = np.sign(data).astype(np.float64)

    else:
        raise ValueError(f"⚠️ '{method}' is not a valid normalization method")

    return normalized_data

def whiten_signal(data: np.ndarray, freq_min: int, freq_max: int, total_channels: int,
                  sampling_frequency: int, axis: int)-> np.ndarray:

    """Performs spectral whitening of all channels. Adapted code from: https://github.com/seismo-live/seismo_live .
    Signals should be adequatly pre-processed.

    Parameters
    ----------
    data : np.ndarray
        Data to spectrally whiten.
    freq_min,freq_max : int, int
        Minimum and maximum of frequency band in which to perform spectral whitening.
    total_channels : int
        Total number of channels.
    sampling_frequency : int
        Sampling frequency of data

    Returns
    -------
    whitened_matrix : np.ndarray
        Whitened data matrix

    """

    whitened_matrix = np.zeros_like(data, dtype="float32")

    for i in tqdm(range(total_channels), desc='Whitening', leave=False):

        channel = data[:, i] if axis == 0 else data[i, :]
        n = len(channel)
        f_range = float(freq_max) - float(freq_min)
        nsmo = int(np.fix(min(0.01, 0.5 * f_range) * float(n) / sampling_frequency))
        f = np.arange(n) * sampling_frequency / (n - 1.0)
        JJ = ((f > float(freq_min)) & (f < float(freq_max))).nonzero()[0]

        # channel FFT
        FFTs = np.fft.fft(channel)
        FFTsW = np.zeros(n, dtype=complex)

        # apodization left
        smo1 = (np.cos(np.linspace(np.pi / 2, np.pi, nsmo + 1)) ** 2)
        FFTsW[JJ[0]:JJ[0] + nsmo + 1] = smo1 * np.exp(1j * np.angle(FFTs[JJ[0]:JJ[0] + nsmo + 1]))

        # boxcar
        FFTsW[JJ[0] + nsmo + 1:JJ[-1] - nsmo] = np.ones(len(JJ) - 2 * (nsmo + 1)) * \
        np.exp(1j * np.angle(FFTs[JJ[0] + nsmo + 1:JJ[-1] - nsmo]))

        # apodization to the right
        smo2 = (np.cos(np.linspace(0.0, np.pi / 2.0, nsmo + 1)) ** 2)
        espo = np.exp(1j * np.angle(FFTs[JJ[-1] - nsmo:JJ[-1] + 1]))
        FFTsW[JJ[-1] - nsmo:JJ[-1] + 1] = smo2 * espo

        # channel IFFT
        whitedata = 2.0 * np.fft.ifft(FFTsW).real

        if axis == 0:
            whitened_matrix[:, i] = np.require(whitedata, dtype='float32')
        else:
            whitened_matrix[i, :] = np.require(whitedata, dtype='float32')

    return whitened_matrix

def signal_spectrum(o_signal: np.ndarray, fs: int, mode: str = "spectrum", pre_processing: bool = True,
    norm: bool = False, order: int = 1, pad: int = 0, nfft: int | None = None,
    nperseg: int | None = None, axis: int = None) -> tuple[np.ndarray, np.ndarray]:

    """Computes the spectrum for a signal, either through FFT or Welch's method.

    Parameters
    ----------
    o_signal : np.ndarray
        Signal to process.
    fs : int
        Sampling frequency of signal.
    mode : str, optional
        ``"spectrum"`` or ``"psd"``.
    pre_processing : bool, optional
        Signal is demeaned, detrended and tapered.
    norm : bool, optional
        If ``True``, result is normalized by maximum value.
    order : int, optional
        Order of polynomial to detrend.
    pad : int, optional
        Number of zeros to add to the signal before and after to increase num of points.
    nfft : int | None, optional
        Number of samples for FFT.
    nperseg : int | None, optional
        Length of each segment for Welch's method.
    axis : int
        Axis along which to compute the spectrum.

    Raises
    ------
    ValueError
        Not a valid mode chosen.

    Returns
    -------
    positive_freqs : np.ndarray
        Frequency vector.
    magnitude : np.ndarray
        Magnitude vector.

    """

    if pre_processing:
        o_signal = filt_preprocess(io_signal=o_signal, order=order, axis=axis)

    if mode == "spectrum":
        if pad > 0:
            if o_signal.ndim == 1:
                o_signal = np.pad(o_signal, (pad-1, pad), mode='constant')
            elif axis == 0:
                o_signal = np.pad(o_signal, ((pad-1, pad),(0, 0)), mode='constant')
            else:
                o_signal = np.pad(o_signal, ((0, 0),(pad-1, pad)), mode='constant')

        n = o_signal.shape[axis] if nfft is None else nfft
        fft = np.fft.fft(o_signal, n=n, axis=axis)
        freq_axis = np.fft.fftfreq(n, 1 / fs)
        positive_freqs = freq_axis[:n//2]

        if o_signal.ndim == 1:
            magnitude = 2/n * np.abs(fft)[:n//2]
        elif axis == 0:
            magnitude = 2/n * np.abs(fft)[:n//2, :]
        else:
            magnitude = 2/n * np.abs(fft)[:, :n//2]

    elif mode == "psd":
        positive_freqs, magnitude = signal.welch(o_signal, fs, nperseg=nperseg, axis=axis)

    else:
        raise ValueError("⚠️ Invalid mode. Choose one of:\n"
                         " - 'spectrum'\n - 'psd'")
    if norm:
        magnitude = magnitude / magnitude.max()

    return positive_freqs, magnitude

def signal_spectrogram(data: np.ndarray, sampling_frequency: int, axis: int,
                       norm: bool)-> tuple[np.ndarray, np.ndarray, np.ndarray]:

	"""Computes spectrogram of signal.

    Parameters
    ----------
    data : np.ndarray
        Signal to process.
    sampling_frequency : int
        Sampling frequency of signal.
    axis : int
        Time axis index.
    norm : bool
        Toggles normalization of results by maximum value.

    Returns
    -------
    f : np.ndarray
        Frequency axis vector.
    t : np.ndarray
        Time axis vector.
    Sxx : np.ndarray
        Spectrogram image.

    """
	nyquist = sampling_frequency/2
	nfft, nperseg = nyquist*2, int(sampling_frequency/5)
	noverlap = int(nperseg/2)
	f, t, Sxx = signal.spectrogram(data, sampling_frequency, nfft=nfft, nperseg=nperseg, noverlap=noverlap)
	Sxx = np.flip(Sxx, axis=axis)
	Sxx = Sxx / Sxx.max(axis=axis) if norm == True else Sxx

	return f, t, Sxx