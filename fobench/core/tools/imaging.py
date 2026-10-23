"""Virtual shot gather and dispersion image computation."""

import pickle as pkl
from pathlib import Path
import numpy as np
from obspy.core import UTCDateTime as UTC
from obspy.signal.cross_correlation import correlate
import pyqtgraph as pg
from pyqtgraph.exporters import ImageExporter
from PyQt5 import QtCore
from fobench.core.tools.signals import taper_signal
from fobench.core.plotting import plotting_pyqt as plot_pyqt

class VSG:
    """Computes a virtual shot gather from a ``Fiber`` instance.

    All channels are cross-correlated with the virtual source channel.
    Positive lags correspond to waves travelling from the virtual source to the
    receiver channel.

    Parameters
    ----------
    fiber : Fiber object
        ``Fiber`` class object, should be pre-processed. Processing is applied
        to ``Fiber`` instance directly.
    max_lag : float
        Maximum lag of the correlation functions in seconds.
    source_channel : int | None
        Channel number used as virtual source. If ``None`` (default),
        middle channel is used.
    temporal_norm : str | None
        Normalization method passed to :meth:`Fiber.normalize`, e.g. ``"1bit"``
        or ``"running mean"``. ``None`` disables normalization.
    ram_window : float
        Window length of the running absolute mean normalization in seconds,
        passed to :meth:`Fiber.normalize`.
    whiten : tuple[float, float] | None
        Frequency band ``(freq_min, freq_max)`` passed to :meth:`Fiber.whiten`.
        ``None`` disables whitening.
    normalize : bool
        Normalizes each correlation function normalization with the overall
        standard deviations. See :func:`obspy.signal.cross_correlation.correlate`.

    Raises
    ------
    ValueError
        ``m̀ax_lag`` longer than record.

    """

    def __init__(self, fiber, max_lag: float, source_channel: int = None,
                 temporal_norm: str = None, ram_window: float = None,
                 whiten: tuple[float, float] = None, normalize: bool = True):

        if source_channel is None:
            source_channel = fiber.channels[fiber.n_channels // 2]
        if max_lag * fiber.sampling_rate >= fiber.n_samples:
            raise ValueError("⚠️ max_lag must be shorter than the record!")

        if whiten is not None:
            fiber.whiten(freq_min=whiten[0], freq_max=whiten[1])
        if temporal_norm is not None:
            fiber.normalize(method=temporal_norm, ram_window=ram_window)

        # metadata
        self.source_channel = int(source_channel)
        self.source_index = int(np.where(fiber.channels == source_channel)[0][0])
        self.channels = np.asarray(fiber.channels).copy()
        self.distances = np.asarray(fiber.distances).copy()
        self.offsets = self.distances - self.distances[self.source_index]
        self.n_channels = fiber.n_channels
        self.sampling_rate = fiber.sampling_rate
        self.dt = fiber.dt
        self.spatial_interval = fiber.spatial_interval
        self.start_time = fiber.start_time
        self.end_time = fiber.end_time
        self.units = fiber.units
        self.fiber_processing = fiber.processing

        self.parameters = {"max_lag": max_lag, "temporal_norm": temporal_norm,
                           "ram_window": ram_window, "whiten": whiten,
                           "normalize": normalize}
        self.created = UTC.utcnow().ctime()

        # correlation
        data = np.moveaxis(fiber.data, (fiber.__axis__("t"), fiber.__axis__("d")), (0, 1))
        data = data.astype(np.float64)
        n_lag = int(round(max_lag * self.sampling_rate))
        source = data[:, self.source_index]
        self.data = np.empty((2 * n_lag + 1, self.n_channels))
        for i in range(self.n_channels):
            self.data[:, i] = correlate(data[:, i], source, n_lag, demean=True,
                                        normalize="naive" if normalize else None,
                                        method="auto")
        self.lags = np.arange(-n_lag, n_lag + 1) * self.dt

    def __str__(self):
        """Defines output of print(VSG)."""
        attributes = ["units", "start_time", "end_time", "source_channel", "n_channels",
                      "spatial_interval", "sampling_rate"]
        values = {attr: getattr(self, attr) for attr in attributes} | self.parameters

        return ("Instance of VSG class\n"
                "virtual shot gather parameters:\n"
                f"{'-' * 65}\n" + '\n'.join(f"{key.ljust(25)} = {value}" for key, value in values.items()))

    __repr__ = __str__

    def plot(self, vmin=None, vmax=None, cmap="seismic", export=None, show=True):
        """Plots the virtual shot gather.
        See :func:`~fobench.core.plotting.plotting_pyqt.plot_2d_distance`.

        Parameters
        ----------
        vmin, vmax : float
            Minimum and maximum limits of colorbar, defaults to the negative and
            positive 95th percentile of the absolute data.
        cmap : str
            Colormap to use.
        export : str
            Filepath to save the plot to.
        show : bool
            Toggle displaying plot.

        Returns
        -------
        plot : pg.PlotItem
            PyQtGraph PlotItem of the VSG.
        win : pg.GraphicsLayoutWidget
            PyQtGraph GraphicsLayoutWidget..

        """
        p95 = np.percentile(np.abs(self.data), 95)
        vmin = -p95 if vmin is None else vmin
        vmax = p95 if vmax is None else vmax

        plot, win = plot_pyqt.plot_2d_distance.__wrapped__( #avoid the decorator
            distances=self.offsets, channels=self.channels - self.source_channel,
            data=self.data.T, y_ticks=self.lags, y_label="Lag [s]",
            x_label="Channel Offset", distance_label="Offset [m]",
            title=f"Virtual Shot Gather (Source: {self.source_channel})",
            cmap=cmap, cbar_label="Correlation Coefficient"
            if self.parameters["normalize"] else "Correlation",
            invert_y=True, vmin=vmin, vmax=vmax)
        plot.addItem(pg.InfiniteLine(pos=0, angle=90,
                                     pen=pg.mkPen("k", width=2, style=QtCore.Qt.DashLine)))
        if export is not None:
            ImageExporter(plot).export(str(export))
        if show:
            win.show()
            pg.exec()
        return plot, win

    def save(self, path: str):
        """Saves the computed ``VSG`` to ``path``."""
        path = Path(path)
        if path.suffix == "":
            path = path.with_suffix(".pkl")
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "wb") as f:
            pkl.dump(self, f)
        return path

    @staticmethod
    def load(path: str):
        """Loads a ``VSG``object from ``path``."""
        with open(path, "rb") as f:
            obj = pkl.load(f)
        if not isinstance(obj, VSG):
            raise TypeError(f"⚠️ {path} is not a VSG object!")
        return obj

    def dispersion(self, c_min: float = 100, c_max: float = 4000, dc: float = 1.0,
                   f_min: float = None, f_max: float = None, normalize: bool = True,
                   taper: float = None, plot_mode: str = "pyqt",
                   export: str = None, show: bool = True, **kwargs):

        """Computes the dispersion image of the virtual shot gather using the
        phase shift method after Park et al. (1998). See :func:`dispersion_image`
        for more details.
        """

        zero = len(self.lags) // 2
        data = (self.data[zero:] + self.data[zero::-1]) / 2

        kwargs.setdefault("title", f"Dispersion Image (Source: {self.source_channel})")
        return dispersion_image(data=data, offsets=np.abs(self.offsets), dt=self.dt,
                                c_min=c_min, c_max=c_max, dc=dc, f_min=f_min, f_max=f_max,
                                axis=0, normalize=normalize, taper=taper, plot_mode=plot_mode,
                                export=export, show=show, **kwargs)


def dispersion_image(data: np.ndarray, offsets: np.ndarray, dt: float, c_min: float = 100,
                     c_max: float = 4000, dc: float = 1.0, f_min: float = None, f_max: float = None,
                     axis: int = 0, normalize: bool = False, taper: float = None,
                     plot_mode: str = "pyqt", export: str = None, show: bool = True, **kwargs):
    """Computes the dispersion image of a shot gather using the phase shift
    method after Park et al. (1998).

    Parameters
    ----------
    data : np.ndarray
        2D array of traces.
    offsets : np.ndarray
        Source-receiver offsets of the traces in m.
    dt : float
        Sampling rate in seconds.
    c_min, c_max : float
        Minimum and maximum phase velocity in m/s.
    dc : float
        Phase velocity increment in m/s.
    f_min, f_max : float | None
        Frequency range in Hz. Defaults to the first non-zero frequency and
        the Nyquist frequency.
    axis : int
        Time axis of ``data``.
    normalize : bool
        Normalizes the dispersion image to the maximum of each frequency if ``True``.
        If ``False`` the image is normalized by the number of traces, i.e. values
        between 0 and 1.
    taper : float | None
        Shape parameter of a Tukey window applied in time, ``None`` disables
        tapering.
    plot_mode : str | None
        ``"pyqt"`` plots the dispersion image with the picked dispersion
        curve, see :func:`~fobench.core.plotting.plotting_pyqt.plot_dispersion`.
        Any other value, e.g. ``None``, disables plotting.
    export : str
        Filepath to save the plot to, only used if plotting is enabled.
    show : bool
        Toggle displaying plot.
    **kwargs
        Arguments passed to :func:`~fobench.core.plotting.plotting_pyqt.plot_dispersion`,
        e.g. ``title``, ``cmap``, ``vmin``, ``vmax``.

    Returns
    -------
    frequencies : np.ndarray
        Frequency vector in Hz.
    velocities : np.ndarray
        Phase velocity vector in m/s.
    image : np.ndarray
        Dispersion image of shape ``(n_velocities, n_frequencies)``.
    curve : np.ndarray
        Phase velocity of the image maximum at each frequency.

    Raises
    ------
    ValueError
        Invalid velocity range.

    """
    data = np.moveaxis(np.asarray(data), axis, 0)
    if not 0 < c_min < c_max:
        raise ValueError("⚠️ velocities must in the right range: 0 < c_min < c_max!")

    data = data - data.mean(axis=0)
    if taper is not None:
        data = taper_signal(data, axis=0, alpha=taper)

    spectrum = np.fft.rfft(data, axis=0)
    frequencies = np.fft.rfftfreq(data.shape[0], dt)
    f_min = frequencies[1] if f_min is None else f_min
    f_max = frequencies[-1] if f_max is None else f_max
    f_mask = (frequencies >= f_min) & (frequencies <= f_max) & (frequencies > 0)
    frequencies, spectrum = frequencies[f_mask], spectrum[f_mask]

    amplitude = np.abs(spectrum)
    spectrum = np.divide(spectrum, amplitude, out=np.zeros_like(spectrum),
                         where=amplitude > 0)

    velocities = np.arange(c_min, c_max + dc/2, dc)
    image = np.empty((len(velocities), len(frequencies)))
    delays = offsets[None, :] / velocities[:, None]
    for i, f in enumerate(frequencies):
        image[:, i] = np.abs(np.exp(2j * np.pi * f * delays) @ spectrum[i])

    if normalize:
        image /= np.where((im_max := image.max(axis=0)) > 0, im_max, 1)
    else:
        image /= len(offsets)
    curve = velocities[np.argmax(image, axis=0)]

    if plot_mode == "pyqt":
        plot_pyqt.plot_dispersion(frequencies, velocities, image, curve,
                                  export=export, show=show, **kwargs)
    return frequencies, velocities, image, curve