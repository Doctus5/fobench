"""Virtual shot gather and dispersion image computation."""

import pickle as pkl
from pathlib import Path
import numpy as np
from obspy.core import UTCDateTime as UTC
from obspy.signal.cross_correlation import correlate
import pyqtgraph as pg
from pyqtgraph.exporters import ImageExporter
from PyQt5 import QtCore
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
        self.fiber_processing = list(fiber.processing)

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
        show : bool, optional
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
            distances=self.distances, channels=self.channels, data=self.data.T,
            y_ticks=self.lags, y_label="Lag [s]",
            title=f"Virtual Shot Gather (Source: {self.source_channel})",
            cmap=cmap, cbar_label="Correlation Coefficient"
            if self.parameters["normalize"] else "Correlation",
            invert_y=True, vmin=vmin, vmax=vmax)
        source_line = pg.InfiniteLine(pos=self.source_channel, angle=90,
                                      pen=pg.mkPen("k", width=2, style=QtCore.Qt.DashLine))
        plot.addItem(source_line)
        x_axis = plot.getAxis("bottom")
        plot.getViewBox().sigXRangeChanged.connect(
            lambda: source_line.setValue(self.source_channel if x_axis.labelText == "Channel"
                                         else self.distances[self.source_index]))
        if export is not None:
            ImageExporter(plot).export(str(export))
        if show:
            win.show()
            pg.exec()
        return plot, win

    def pickle(self, path: str):
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


def dispersion_image():
    pass