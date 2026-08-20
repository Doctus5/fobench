"""Contains all functionality related to plotting using matplotlib, i.e.
whenever plot_mode is set to 'mpl'."""

import numpy as np
import datetime as datetime
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.dates import num2date

"""Helper Functions"""

class PrecisionDateFormatter(ticker.Formatter):
    """Extends the `matplotlib.ticker.Formatter` class to allow for millisecond
    precision when formatting a tick (in days since the epoch) with a
    `~datetime.datetime.strftime` format string. Adapted from StackOverflow.
    """

    def __init__(self, fmt: str, precision: int = 3, tz: str = None):
        """
        Parameters
        ----------
        fmt : str
            `~datetime.datetime.strftime` format string.
        precision : int
            Necessary precision.
        tz : str
            Timezone info.
        """

        self.num2date = num2date
        self.fmt = fmt
        self.tz = tz if tz is not None else datetime.timezone.utc
        self.precision = precision

    def __call__(self, x, pos=0):
        if x == 0:
            raise ValueError("DateFormatter found a value of x=0, which is "
                             "an illegal date; this usually occurs because "
                             "you have not informed the axis that it is "
                             "plotting dates, e.g., with ax.xaxis_date()")

        dt = self.num2date(x, self.tz)
        ms = dt.strftime("%f")[:self.precision]

        return dt.strftime(self.fmt).format(ms=ms)


"""matplotlib plotting Functions"""

def gen_DAS_plot(data: np.ndarray = None, t: np.ndarray = None, channels: list = None,
                 units_y: str = None, figsize: list[float, float] | tuple[float, float] = None,
                 title: str = None, cmap:str = "seismic", file_name: str = None,
                 add_data: np.ndarray = None, show: bool = True, **kwargs)-> plt.Figure:

    """Generic FOS data plot.

    Parameters
    ----------
    data : np.ndarray
        Data to plot.
    t : np.ndarray
        Timestamps in matplotlib format.
    channels : list[int | float]
        List of channel numbers or distance values.
    units_y : str
        The units of the amplitude values.
    figsize : list[float, float] | tuple[float, float]
        Controls figure size
    title : str
        Title of plot
    cmap : str, optional
        matplotlib colormap to use.
    file_name : str
        If not ``None``, plot will be saved at given location.
    add_data : np.ndarray
        Additional data to plot on top of FOS data.
    show : bool
        Toggles display of figure.
    **kwargs:
        Additional arguments passed to ``imshow``

    Returns
    -------
    fig : plt.Figure
        Figure instance.
    """

    fig, ax = plt.subplots(figsize=figsize)
    fig.autofmt_xdate()
    extent=(channels[0], channels[-1]+1,t[-1],t[0])
    cm = ax.imshow(data, cmap=cmap, extent=extent,
                   aspect="auto", interpolation="none", **kwargs)
    plt.colorbar(cm, label=units_y.title())

    if add_data is not None:
        ax.imshow(add_data, cmap="binary", extent=extent, aspect="auto", alpha=add_data)

    ax.yaxis_date()
    precision = str(datetime.timedelta(days=(t[1]-t[0])).total_seconds())[::-1].find(".")
    ax.yaxis.set_major_formatter(PrecisionDateFormatter("%H:%M:%S.{ms}", precision))
    ax.set_ylabel("Time", fontsize=15)
    ax.set_xlabel("Channel", fontsize=15)
    ax.set_title(title, fontsize=20)
    ax.tick_params(axis="x", bottom=True, top=True, labelbottom=True,
                   labelleft=True, labeltop=True, rotation=0)

    if file_name is not None:
        fig.savefig(file_name, transparent=False, bbox_inches="tight", pad_inches=0)
    if show:
        plt.show()

    return fig


def simple_plot(data: np.ndarray, t: np.ndarray, channel: list[str | int | float] = [""],
                units_y: str = None, max_value: float = None,
                figsize: list[float, float] | tuple[float, float] = None, title: str = None,
                file_name: str = None, show: bool = True, **kwargs) -> plt.Figure:

    """Plots a simple timeseries.

    Parameters
    ----------
    data : np.ndarray
        Timeseries to plot.
    t : np.ndarray
        Timestamps in matplotlib format.
    channel : list[str | int | float]
        List of channel numbers or description.
    units_y : str
        The units of the amplitude values.
    max_value : float
        Limits the maximum value of y-axis, if not ``None``, limits will be set to
        -``max_value`` to ``max_value``
    figsize : list[float, float] | tuple[float, float]
        Controls figure size
    file_name : str
        If not ``None``, plot will be saved at given location.
    title : str
        Title of plot
    show : bool
        Toggles display of figure.
    **kwargs:
        Additional arguments passed to ``plot``

    Returns
    -------
    fig : plt.Figure
        Figure instance.
    """

    fig, ax = plt.subplots(1, 1, sharex=True, gridspec_kw={"hspace": 0.3},
                           figsize=figsize)
    fig.autofmt_xdate()

    for ch_data, label in zip(data, channel):
        ax.plot(t, ch_data, linewidth=0.7, label=label, **kwargs)
    min_val, max_val = (-max_value, max_value) if max_value is not None else (data.min(), data.max())
    ax.set_ylim(min_val * 1.2, max_val * 1.2)

    ax.xaxis_date()
    precision = str(datetime.timedelta(days=(t[1]-t[0])).total_seconds())[::-1].find(".")
    ax.xaxis.set_major_formatter(PrecisionDateFormatter("%H:%M:%S.{ms}", precision))
    ax.set_xlim(t[0],t[-1])
    ax.legend(loc=1)
    ax.set_ylabel(units_y.title(), fontsize=15)
    ax.set_xlabel("Time", fontsize=15)
    ax.set_title(title, fontsize=20)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-2,2))

    if file_name is not None:
        fig.savefig(file_name, transparent=False, bbox_inches="tight", pad_inches=0)
    if show:
        plt.show()

    return fig


def mpl_fx_plot(spec_matrix: np.ndarray = None, freqs: np.ndarray = None,
                x: list = None, units_y: str = None,
                figsize: list[float, float] | tuple[float, float] = None,
                cmap:str = "viridis", title: str = None, file_name: str = None,
                show: bool = True, **kwargs) -> None:

    """Plots the fx-plot, i.e. frequency content over distance.

    Parameters
    ----------
    spec_matrix : np.ndarray
        fx data.
    freqs : np.ndarray
        Frequency axis values.
    x : list
        Channel numbers or distance values.
    units_y : str
        Amplitude units.
    figsize : list[float, float] | tuple[float, float]
        Controls figure size.
    cmap : str, optional
        matplotlib colormap to use for the spectrogram.
    title : str, optional
        Title of plot.
    file_name : str, optional
        If not ``None``, plot will be saved at given location.
    show : bool
        Toggles display of figure.
    **kwargs
        Addtional arguments passed to ``imshow``, e.g. ``vmin``, ``vmax``.

    Returns
    -------
    fig : plt.Figure
        Figure instance.

    """

    fig, ax = plt.subplots(figsize=figsize)

    ax.set_ylabel("Frequency [Hz]", fontsize=15)
    ax.set_xlabel("Channel", fontsize=15)
    ax.set_title(title, fontsize=20)

    extent = (x[0], x[-1]+1, freqs[0], freqs[-1])
    cm = ax.imshow(spec_matrix, cmap=cmap, extent=extent, aspect="auto",
                   interpolation="none", **kwargs)
    plt.colorbar(cm, ax=ax, label=units_y)
    ax.tick_params(axis="x", bottom=True, top=True, labelbottom=True, labelleft=True, labeltop=True, rotation=0)
    if file_name is not None:
        fig.savefig(file_name, transparent=False, bbox_inches="tight", pad_inches=0)
    if show:
        plt.show()

    return fig


def simple_spectrogram(data: np.ndarray = None, freq: np.ndarray = None,
                       t: np.ndarray = None, units_y: str = None,
                       figsize: list[float, float] | tuple[float, float] = None,
                       trace: bool = False, cmap: str = "viridis", title: str = "",
                       file_name: str = None, show: bool = True,
                       freq_lim: list[float, float] | tuple[float, float] = None,
                       **kwargs) -> plt.Figure:

    """Plots a simple spectrogram, optionally together with the time domain data.

    Parameters
    ----------
    data : np.ndarray
        Spectrogram data.
    freq : np.ndarray
        Frequency axis values.
    t : np.ndarray
        Timestamps in matplotlib format.
    units_y : str
        Amplitude units.
    figsize : list[float, float] | tuple[float, float]
        Controls figure size.
    trace : bool
        Toggles plotting of the time domain data in a second panel.
    cmap : str, optional
        matplotlib colormap to use for the spectrogram.
    title : str, optional
        Title of plot.
    file_name : str, optional
        If not ``None``, plot will be saved at given location.
    freq_lim : list[float, float] | tuple[float, float], optional
        Limits for the displayed frequency range in spectrogram.
    show : bool
        Toggles display of figure.
    **kwargs
        Addtional arguments passed to ``imshow``, e.g. ``vmin``, ``vmax``.

    Returns
    -------
    fig : plt.Figure
        Figure instance.

    """

    extent = (t[0], t[-1], freq[0], freq[-1])
    precision = str(datetime.timedelta(days=(t[1]-t[0])).total_seconds())[::-1].find(".")

    if trace is not None:
        fig, ax = plt.subplots(2,1, sharex=True,
                               gridspec_kw={"hspace": 0.2,"height_ratios":[3,1]},
                               figsize=figsize)
        cm = ax[0].imshow(data, cmap=cmap, extent=extent, aspect="auto",
                          interpolation="none", **kwargs)
        plt.colorbar(cm, ax=ax[0], label=units_y.title(), orientation="horizontal",
                     aspect=40, pad=0.02)
        ax[0].xaxis.set_major_formatter(PrecisionDateFormatter("%H:%M:%S.{ms}", precision))
        ax[0].set_xlim(t[0],t[-1])
        ax[0].set_title(title, fontsize=20)
        ax[0].set_ylabel("Frequency [Hz]", fontsize=10)

        ax[1].plot(t, trace, c="black", linewidth=0.7)
        ax[1].set_ylabel(units_y.title(), fontsize=10)
        ax[1].set_xlabel("Time", fontsize=15)

        if freq_lim is not None:
            ax[0].set_ylim(freq_lim[0],freq_lim[1])

        ax[1].tick_params(axis="x", bottom=True, top=False, labelbottom=True,
                          labelleft=True, labeltop=False, labelrotation=25)

    else:
        fig, ax = plt.subplots(figsize=figsize)
        fig.autofmt_xdate()

        cm = ax.imshow(data, cmap=cmap, extent=extent, aspect="auto",
                       interpolation="none", **kwargs)
        plt.colorbar(cm, ax=ax, label=units_y.title())
        ax.xaxis.set_major_formatter(PrecisionDateFormatter("%H:%M:%S.{ms}", precision))
        ax.set_xlim(t[0],t[-1])
        ax.set_title(title, fontsize=20)
        ax.set_ylabel("Frequency [Hz]", fontsize=15)
        ax.set_xlabel("Time", fontsize=15)
        ax.tick_params(axis="x", bottom=True, top=False, labelbottom=True,
                       labelleft=True, labeltop=False, labelrotation=25)
        if freq_lim is not None:
            ax.set_ylim(freq_lim[0],freq_lim[1])

    if file_name is not None:
        fig.savefig(file_name, transparent=False, bbox_inches="tight", pad_inches=0)
    if show:
        plt.show()

    return fig


def simple_spectrum(spectra: np.ndarray = None, freqs: np.ndarray = None,
                    channels: list = None, y_units: str = None, legend: bool = True,
                    figsize: tuple[float, float] = None, file_name: str = None,
                    title: str = None, show: bool = True, **kwargs) -> plt.Figure:

    """Plots a simple spectrum for one or multiple channels.

    Parameters
    ----------
    spectra : np.ndarray
        Spectrum data to plot.
    freqs : np.ndarray
        Frequency axis values.
    channels : np.ndarray
        Channel numbers.
    y_units : str
        Units of the spectral amplitude.
    legend : bool
        Toggles legend in plot.
    figsize : tuple[float, float]
        Controls figure size.
    file_name : str
        If not ``None``, plot will be saved at given location.
    title : str
        Title of plot.
    show : bool
        Toggles display of figure.
    **kwargs
        Additional arguments passed to ``plot``.

    Returns
    -------
    fig : plt.Figure
        Figure instance.

    """

    fig, ax = plt.subplots(figsize=figsize)

    ax.ticklabel_format(axis="y", style="sci", scilimits=(-2,2))
    ax.set_xlabel("Frequency [Hz]", fontsize=15)
    ax.set_ylabel(y_units.title(), fontsize=15)
    ax.set_title(title, fontsize=20)
    ax.set(xlim=(freqs[0], freqs[-1]), ylim=(0, spectra.max() * 1.1))
    ax.grid()

    for channel, spectrum in zip(channels, spectra):
        ax.plot(freqs, spectrum, label=channel, **kwargs)

    if legend:
        ax.legend()
    if file_name is not None:
        fig.savefig(file_name, transparent=False, bbox_inches="tight", pad_inches=0)
    if show:
        plt.show()

    return fig


def plot_record_section(signals: np.ndarray, t: np.ndarray, channels: np.ndarray,
                        date: str, show: bool = True) -> plt.Figure:

    """Plots a record section for selected channels.

    Parameters
    ----------
    signals : np.ndarray
        Channel data to plot.
    t : np.ndarray
        Timestamps in matplotlib format.
    channels : np.ndarray
        Channel numbers.
    date : str
        The start date of data.
    show : bool
        Toggles display of figure.

    Returns
    -------
    fig : plt.Figure
        Figure instance.

    """

    num_stations = len(channels)
    scaling_factor = 1.5 / np.max(np.abs(signals))

    fig, ax = plt.subplots(figsize=(10, 8))

    for i in range(num_stations):
        ax.plot(t, signals[:, i] * scaling_factor +
                i, color="black", linewidth=0.7)

    ax.set_yticks(np.arange(num_stations))
    ax.set_yticklabels(channels)
    ax.invert_yaxis()

    ax.xaxis_date()
    precision = str(datetime.timedelta(
        days=(t[1]-t[0])).total_seconds())[::-1].find(".")
    ax.xaxis.set_major_formatter(
        PrecisionDateFormatter("%H:%M:%S.{ms}", precision))
    ax.set_xlim(t[0], t[-1])
    ax.grid(color="gray", linestyle="--", alpha=0.8)
    ax.set(xlabel="Time [s]", ylabel="Channel")
    ax.set_title(f"Record Section for {date}", fontsize=20)
    if show:
        plt.show()

    return fig


def plot_acfs(acfs: np.ndarray, distances: list | np.ndarray, fs: int,
              max_shift: int, **imshow_kwargs) -> None:

    """Plots the autocorrelation profile

    Parameters
    ----------
    acfs : np.ndarray
        Array containing the autocorrelations.
    distances : list | np.ndarray
        Channel distances.
    fs : int
        Sampling frequency of data.
    max_shift : int
        Maximum shift chosen.
    **imshow_kwargs :
        Addtional arguments passed to ``imshow``, e.g. ``vmin``, ``vmax``.

    Returns
    -------
    None

    """

    extent = [distances[0], distances[-1], max_shift/fs, 0]
    y_label, x_label = "Lag/TWT [s]", "Optical Distance [m]"
    im = plt.imshow(acfs, aspect="auto", origin="upper", extent=extent,
                    **imshow_kwargs)
    vmin = imshow_kwargs.get("vmin")
    vmax = imshow_kwargs.get("vmax")

    if vmin is not None and vmax is not None:
        extend = "both"
    elif vmin is not None:
        extend = "min"
    elif vmax is not None:
        extend = "max"
    else:
        extend = "neither"

    plt.colorbar(im, label="Auto-correlation Coefficient", extend=extend)
    plt.xlabel(x_label, fontsize=16)
    plt.ylabel(y_label, fontsize=16)
    plt.gca().tick_params(axis="both", labelsize=16)
    plt.tight_layout()