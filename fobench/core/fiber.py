"""``Fiber`` class for reading, manipulating and visualizing fiber optic sensing data.

:Authors:
    - Sergio Diaz-Meza
    - Jonas Pätzel

:Contributors:
    - Christopher Wollin

"""

import sys
import copy
from warnings import warn

import numpy as np
from obspy.core import UTCDateTime as UTC

import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets

from fobench.core.tools import file_io, utils, filters, signals, wavefield
from fobench.core.plotting import plotting_mpl as plot
from fobench.core.plotting import plotting_pyqt as plot_pyqt
from fobench.core.plotting.pyqt_viewer import Viewer


class Fiber(object):
    """The base class of FoBench, reading in data and metadata.

    Note
    ----
    - Most of the methods perform changes within the class permanently.
      It can be useful to make a copy of the Fiber instance with the
      ``.copy()`` method before performing any processing or changes.
    - For most methods the ``plot_mode`` parameter can either be
      ``"pyqt"`` or ``"mpl"``; to not generate a plot, set it to anything
      else, e.g. ``None``.
    """

    def __init__(self, filepath: str, company: str = "", range_ch: tuple[int,int] = None, sensing: str = "das",
                 load_data: bool = True, show_progress: bool = True, storage_opts = None):

        """

        Parameters
        ----------
        filepath : str
            Path to file to read.
        company : str
            Interrogator manufacturer. One of ``"silixa"``, ``"febus"``, ``"aragon"``,
            ``"quantx"``, ``"asn"``, ``"terra15"`` ``bam`` and ``"sintela"``
        range_ch : tuple(int, int)
            Range of channels to load.
        sensing : str
            Fiber optic sensing technology. Currently only ``"das"``.
        load_data : bool
            If ``False``, only metadata is loaded.
        show_progress : True
            Show progress bar when loading data.
        storage_opts :
            -

        """

        if not company:
            raise ValueError(
                "\nNo company provided! Please choose one of:\n"
                " -'silixa'\n -'febus'\n -'bam'\n -'aragon'\n -'quantx'\n -'asn'\n"
                " -'terra15'\n -'sintela'"
            )

        self.__filepath__ = [filepath]
        self.__storage_opts__ = storage_opts

        self.company = company.lower()
        self.format = filepath.split(".")[-1]

        self.attributes = file_io.read_data(self.__filepath__[0], self.company, range_ch,
                        self.format, load_data=load_data, show_progress=show_progress, storage_opts=storage_opts)

        self.__basefile__ = self.attributes["basefile"] # changed to the structure of the file
        self.fiber = self.attributes["fiber"]
        self.properties = self.attributes["properties"] # all metadata of input file
        self.channels = self.attributes["chans"] # list of channels as array
        self.total_channels = self.attributes["total_channels"]
        self.sampling_frequency = self.attributes["sampling_frequency"] # sampling rate of the data.
        self.o_sampling_frequency = self.attributes["o_sampling_frequency"] if self.attributes["o_sampling_frequency"] != None else self.attributes["sampling_frequency"] # original sampling frequency. Important for conversion factor.
        self.dt = 1 / self.sampling_frequency # calculated time step.
        self.start_time = self.attributes["start_time"] # start time of the data in file.
        self.end_time = self.attributes["end_time"] # end time of the data in file.
        self.spatial_interval = self.attributes["spatial_interval"] # channel spacing or spatial interval between channels [m].
        self.time_length = self.end_time - self.start_time
        self.num_points = self.attributes["num_points"] # int(self.time_length/self.dt)
        self.gauge_length = self.attributes["gauge_length"] # gauge length used in the measurement [m].
        self.channel_offset = self.attributes["channel_offset"] # offset where measurement started. It will not always record at channel 0 or distance 0.
        self.data = self.attributes["data"]
        self.corrected = False
        self.sensing = sensing
        self.units = self.attributes["units"]
        self.conv_factor = self.attributes["conv_factor"] # Extra variables (ONLY FOR ASN HDF5)
        self.processing = [{"instance creation" : UTC.utcnow().ctime()}]
        self.distances = (self.channels - self.channel_offset) * self.spatial_interval

        self.ch_coord = None # coordinates of channels, requires more input ot be filled

    """Internat Methods"""

    def __str__(self):
        """Defines output of print(Fiber)."""
        attributes = ["units", "start_time", "end_time", "num_points", "total_channels",
                    "spatial_interval", "sampling_frequency", "gauge_length"]

        return ("Instance of Fiber class\n"
                "recording parameters:\n"
                f"{'-' * 65}\n"+ '\n'.join(f"{attr.ljust(25)} = {getattr(self, attr)}" for attr in attributes))

    __repr__ = __str__

    def __iadd__(self, other):
        """Defines output of Fiber + Fiber."""
        if not isinstance(other, Fiber):
            raise TypeError("Object to add must be instance of Fiber class")

        return self.concatenate(other, fill_gaps=0)

    def __axis__(self, dim):
        """Mapping of axis order, dim can be 't' (time) or 'd' (distance)."""

        return {"t": 0, "d": 1}[dim]

    """Standard methods"""

    def metadata(self, meta_dict=False):
        """Prints out metadata, optionally returns all metadata as dictionary."""

        if meta_dict:
            metainfo = {key: value for key, value in vars(self).items() if not key.startswith("__")}
            return metainfo

        for prop, value in self.properties.items():
                print(f"{prop} = {value}")

    def copy(self):
        """Returns a deep copy the Fiber class object in its current state."""
        self.__dict__.pop("_viewer", None)

        return copy.deepcopy(self)

    def instr_correct(self, target="strain-rate", terra15_gl=None):
        """Performs instrument correction and data conversion for various instrument types.
        See :func:`~fobench.core.utils.instr_corr`."""
        if not self.corrected:
            (self.data, self.units, self.channels,
            self.total_channels, self.gauge_length, self.distances) = utils.instr_corr(self.data, vars(self),
                                    target=target, terra15_gl=terra15_gl, axis=self.__axis__("d"))
            # self.distances = [(num + self.channel_offset) * self.spatial_interval for num in self.channels_num]
            self.corrected = True
            return self

        warn("Instrument correction has already been applied, doing nothing...")
        return self

    def trim(self, t0=None, tf=None):
        """Cuts data between given start and end times, t0 and tf can be ``UTCDateTime``
        or ISOformat style strings.
        """
        t_axis = self.__axis__("t")
        data, start_time, end_time = utils.trim_time(t0=t0, tf=tf, data=self.data,
                                                     times=self.times(), start_time=self.start_time,
                                                     end_time=self.end_time, axis=t_axis)

        self.data, self.start_time, self.end_time, self.time_length, self.num_points = (data,
                start_time, end_time, end_time-start_time, data.shape[t_axis])

        return self

    def restrict_channels(self, ch0, chf):
        """Trims data in space, between ch0 and chf, a single channel is returned
        when ch0 = chf, updates all class attributes.
        """
        d_axis = self.__axis__("d")
        ch0, chf = int(min(ch0, chf)), int(max(ch0, chf))
        channels_list = self.channels.tolist()
        ch0, chf = channels_list.index(ch0), channels_list.index(chf)
        self.data = self.data[:,ch0:chf+1] if d_axis == 1 else self.data[ch0:chf+1,:]
        self.channels = self.channels[ch0:chf+1]
        self.distances = self.distances[ch0:chf+1]
        self.total_channels = len(self.channels)

        return self

    def append_coord(self, n_ch, x_ch, y_ch, z_ch):
        """Attaches channel coordinates for later plotting. Takes 1D arrays of
        channel number (n_ch), longitude and latitude (x_ch and y_ch) and elevation in m (z_ch).
        """
        coords = [n_ch,
                  np.zeros_like(n_ch) if x_ch is None else x_ch,
                  np.zeros_like(n_ch) if y_ch is None else y_ch,
                  np.zeros_like(n_ch) if z_ch is None else z_ch]

        self.ch_coord = np.column_stack(coords)

        return self

    def georeference(self, n_ch, x_ch, y_ch, z_ch, system="decimal", err=None):
        """Takes known channel locations, e.g. from tap tests and interpolates channel locations
        inbetween, attaches new coordinates.
        takes 1D arrays of channel number (n_ch), longitude and latitude (x_ch and y_ch)
        and elevation in m (z_ch), coordinate system can be for lon and lat can be "decimal" or "utm"
        "err" is maximum accepted interpolation error between original metadata location and new interpolated
        location
        """
        x_ch = np.zeros(n_ch.size) if x_ch is None else x_ch
        y_ch = np.zeros(n_ch.size) if y_ch is None else y_ch
        z_ch = np.zeros(n_ch.size) if z_ch is None else z_ch

        n_ch, x_ch, y_ch, z_ch = utils.interpolate_channels(n_ch, x_ch, y_ch,
                                                      z_ch, system, err, self.spatial_interval)
        self.append_coord(n_ch, x_ch, y_ch, z_ch)

        return self

    def get_data(self, channel=None):
        """Returns data, similar to Fiber.data but with option to return only a specified channel.
        """
        if channel is not None:
            index = self.channels.tolist().index(int(channel))
            d_axis = self.__axis__("d")
            return self.data[:,index] if d_axis == 1 else self.data[index,:]

        return self.data

    def times(self, time_type="UTCDateTime"):
        """Returns array of sample times, can be ``'UTCDateTime'``, ``'isoformat'``,
        ``'datetime64'``, ``'matplotlib'`` or ``"unix"``. See :func:`~fobench.core.tools.utils.return_times`.
        """
        return utils.return_times(self, time_type)

    def concatenate(self, input_das=None, fill_gaps=0):
        """Concatenates two ``Fiber`` class objects assuming they have the same
        sampling rate. Channel order is determined automatically, gaps are filled
        with ``"fill_gaps"`` value.
        """

        axis = self.__axis__("t")
        first, second = (self, input_das) if self.start_time <= input_das.start_time else (input_das, self)
        self.data = np.concatenate((first.data, second.data), axis=axis)
        self.start_time, self.end_time = first.start_time, second.end_time
        self.time_length = self.end_time - self.start_time
        self.num_points = self.data.shape[axis]
        self.__filepath__.extend(input_das.__filepath__)

        return self

    def to_traces(self, t_type="obspy"):
        """Returns channels as ``'obspy'`` or ``'pyrocko'`` streams. See
        :func:`~fobench.core.tools.utils.to_traces`.
        """

        return utils.to_traces(self, t_type)

    def to_xarray(self, name=None, use_distance=True):
        """Returns the Fiber class as xarray.DataArray object. See
        :func:`~fobench.core.tools.utils.to_xarray`.
        """

        return utils.to_xarray(self, name, use_distance)


    def write(self, save_path=""):
        """Save data of Fiber in a new file in its the original format. See
        :func:`~fobench.core.tools.file_io.write_data`.
        """
        if isinstance(self.__basefile__, str):
            self.__basefile__ = file_io.scan_template(self.__basefile__,
                                                      company=self.company,
                                                      format=self.format,
                                                      storage_opts=self.__storage_opts__)

        file_io.write_data(self, filepath=save_path, company=self.company)

    """Signal Processing methods"""

    @utils._update_processing
    def spatial_resample(self, rs_type=None):
        """Modifies spatial sampling of the data by adding or removing channels.
        ``"upsampling"`` adds a channel between each channel pair by interpolating the values.
        ``"downsampling"`` removes every second channel.
        See :func:`~fobench.core.tools.utils.spatial_upsampling` and .:func:`~fobench.core.tools.utils.spatial_downsampling`
        """
        if rs_type in ["upsampling", "upsample"]:
            self.data, self.channels = utils.spatial_upsampling(self)
            self.spatial_interval /= 2
        elif rs_type in ["downsampling", "downsample"]:
            self.data, self.channels = utils.spatial_downsampling(self)
            self.spatial_interval *= 2
        else:
            raise ValueError(f"\nInvalid resample type: '{rs_type}'. Choose on of:\n"
                    " -'upsampling'\n -'downsampling'")
        self.total_channels = len(self.channels)

        return self

    @utils._update_processing
    def detrend(self, order=1, dim="t"):
        """Detrend signal with specified order polynomial.See :func:`~fobench.core.tools.signals.detrend_signal`.
        """
        axis = self.__axis__(dim)
        self.data = signals.detrend_signal(self.data, order=order, axis=axis)

        return self

    @utils._update_processing
    def demean(self, dim="t"):
        """Remove mean of signal along specified dimension.See :func:`~fiber.core.tools.signals.demean_signal`.
        """
        axis = self.__axis__(dim)
        self.data = signals.demean_signal(self.data, axis=axis)

        return self

    @utils._update_processing
    def taper(self, alpha=0.05, dim="t", detaper=False):
        """Tapers data see :func:`~fobench.core.tools.signals.taper_signal.`
        """
        axis = self.__axis__(dim)
        self.data = signals.taper_signal(data=self.data, axis=axis,
                                   alpha=alpha, detaper=detaper)

        return self

    @utils._update_processing
    def decimate(self, new_freq=None, dim="t", f_type="fir-remez"):
        """Decimates data to new sampling frequeny or spatial interval (with prefilter).
        Target frequency should divide original sampling frequency evenly. Filter options
        are ``"fir-remez"`` or ``None`` for SciPys default anti-aliasing order 8
        Chebyshev Type I filter. See :func:`~fobench.core.tools.filters.decimate`.

        Warning
        -------
        Careful when decimating using factors >= 13, it is preferable to call
        decimation twice instead! (See :func:`~scipy.signal.decimate`)
        """
        axis = self.__axis__(dim)

        if dim == "t":
            if new_freq is None:
                raise ValueError("new_freq (temporal) must be provided as a positive number in Hz")
            if new_freq <= 0:
                raise ValueError(f"new_freq must be > 0 Hz, got {new_freq}")
            if new_freq > self.sampling_frequency:
                raise ValueError(f"new_freq ({new_freq} Hz) cannot exceed current sampling frequency ({self.sampling_frequency} Hz)")

            if self.sampling_frequency % new_freq != 0:
                warn(f"Decimation to {new_freq} Hz not possible! Decimating to {self.sampling_frequency / int(self.sampling_frequency / new_freq)} Hz instead")
            down_factor = int(self.sampling_frequency / new_freq)
            new_freq = self.sampling_frequency / down_factor

            self.data = filters.decimate(data=self.data, factor=down_factor, f_type=f_type, axis=axis)

            self.sampling_frequency  = new_freq
            self.dt = 1 / self.sampling_frequency
            self.num_points = self.data.shape[axis]

        elif dim == "d": # filtering and decimate spatially
            new_dx = new_freq
            if new_dx is None:
                raise ValueError("new_freq (spatial) must be provided as a positive number in meters")
            if new_dx <= self.spatial_interval:
                raise ValueError(f"target spatial interval must be larger than {self.spatial_interval} meters")

            factor = int(new_freq / self.spatial_interval)
            if not np.isclose(self.spatial_interval * factor, new_dx):
                warn(f"Spatial interval of {new_dx} not possible. Therefore, the new spatial interval would be of {self.spatial_interval*factor} meters")

            self.data = filters.decimate(data=self.data, factor=factor, f_type=f_type, axis=axis)
            self.spatial_interval *=  factor
            self.channels = self.channels[::factor]
            self.distances = (self.channels -self.channel_offset) * self.spatial_interval
            self.total_channels = len(self.channels)

        return self

    @utils._update_processing
    def normalize(self, method="absolute max", dim="t", ram_window=None):
        """Normalize data. Methods are ``"absolute max"``, ``"trace max"``,
        ``"running mean"`` and ``"1bit"``.
        See :func:`~fiber.core.tools.signals.normalize_signal`.
        """
        axis = self.__axis__(dim)
        self.data = signals.normalize_signal(self.data, method=method,
                                       ram_window=ram_window, axis=axis,
                                       fs=self.sampling_frequency, total_channels=self.total_channels,
                                       num_points=self.num_points)

        return self

    @utils._update_processing
    def whiten(self, freq_min=0.01, freq_max=100, dim="t"):
        """Spectral whitening of data.
        See :func:`~fobench.core.tools.signals.whiten_signal`.
        """
        axis = self.__axis__(dim)
        if not any("filter" in preprocessing for preprocessing in self.processing):
                      warn("Data has possibly not been filtered before whitening! Check"
                     "preprocessing and results carefully!\ncontinuing...")
        self.data = signals.whiten_signal(data = self.data, freq_min=freq_min, freq_max=freq_max,
                                    sampling_frequency=self.sampling_frequency,
                                    total_channels=self.total_channels, axis=axis)

        return self

    @utils._update_processing
    def filter(self, f_type=None, freq=None, pre_process=True, alpha=0.05, order=1, sym=True, dim="t",
            **options):
        """Filters data using specified filter, based on Obspy.signal.filter module.
        For frequency filters, data is optionally pre-processed.
        For ``'bandpass'`` and ``'bandstop'`` options, ``freq`` must be tuple.
        See :func:`~fobench.core.tools.filters.point_filter` for all filter options.
        """

        axis = self.__axis__(dim)
        if pre_process and f_type != "median":
            self.preprocess(alpha=alpha, sym=sym, order=order, dim=dim)
        df = self.sampling_frequency if dim == "t" else 1/self.spatial_interval
        self.data = filters.point_filter(f_type=f_type, data=self.data,
                                  df=df, freq=freq, axis=axis, **options)

        return self

    @utils._update_processing
    def preprocess(self, alpha=0.05, order=1, sym=True, dim="t", steps=(True, True, True)):
        """Performs demeaning, detrending and tapering.
        See :func:`~fobench.core.tools.signals.filt_preprocess`.
        """
        axis = self.__axis__(dim)
        self.data = signals.filt_preprocess(io_signal=self.data, order=order,
                                      alpha=alpha, sym=sym, axis=axis, steps=steps)

        return self

    @utils._update_processing
    def fk_filter(self, bands=[{}], propagation="both", alpha=0.3, plot_mode=None,
                  verbose=False, results=False, mode="pass"):
        """Applies frequency wavenumber filter to data.
        See :func:`~fobench.core.tools.filters.fk_filter`.
        """
        out = filters.fk_filter(data=self.data, dt=self.dt, dx=self.spatial_interval,
                                bands=bands, propagation=propagation, alpha=alpha,
                                plot_mode=plot_mode, verbose=verbose, mode=mode,
                                t_axis=self.__axis__("t"), d_axis=self.__axis__("d"))
        self.data = out[0] if verbose else out
        if results:
            return (out[0], out[1], out[2]) if verbose else (out[0])

        return self

    @utils._update_processing
    def integrate(self, dim="t", taper=True):
        """Integrates data using cum-trapezoids method, tapers data by default.
        See :func:`~fobench.core.tools.signals.integrate_signal`.
        """
        axis = self.__axis__(dim)
        dx = self.dt if dim == "t" else self.spatial_interval
        if taper:
            self.taper(dim=dim)
        self.data = signals.integrate_signal(data=self.data, dx=dx, axis=axis)

        return self

    @utils._update_processing
    def differentiate(self, method="gradient", dim="t"):
        """Differentiates data in space or time using ``'gradient'`` or ``'diff'`` method.
        See :func:`~fobench.core.tools.signals.differentiate_signal`.
        """
        if method not in ["gradient", "diff"]:
            raise ValueError(f"\nInvalid method: '{method}'. Choose on of:\n"
                    " -'gradient'\n -'diff'")
        axis = self.__axis__(dim)
        dx = self.dt if dim == "t" else self.spatial_interval
        self.data = signals.differentiate_signal(self.data, method=method,
                                           axis=axis, dt=dx)

        return self

    def SNR(self, dim="t", results=False, plot_mode="pyqt"):
        """Computes signal to noise ratio, defined as ratio between mean and standard
        deviation of signal.
        """
        axis = self.__axis__(dim)
        snr = self.data.mean(axis=axis) / self.data.std(axis=axis)
        if plot_mode == "mpl":
            warn("⚠️ matplotlib plotting not implemented for this method, "
                 "plotting using pyqtgraph instead")
            plot_mode = "pyqt"

        if plot_mode == "pyqt":
            plot_pyqt.plot_distance(distances=self.distances, channels_num=self.channels,
                                       data=snr, y_label="SNR [-]", title="SNR Profile")
        if results:
            return snr

    def rmsa(self, window=None, dim="t", plot_mode="pyqt", results=False,
          vmin=None, vmax=None):
        """Computes root mean square amplitude for record, dependign on dimension,
        window is either seconds ("t") or number of channels ("d").
        See :func:`~fobench.core.tools.wavefield.rmsa`.
        """
        axis = self.__axis__(dim)
        if window is not None and dim == "t": window =  window*self.sampling_frequency
        if plot_mode == "mpl":
            warn("⚠️ matplotlib plotting not implemented for this method, "
                 "plotting using pyqtgraph instead")
            plot_mode = "pyqt"
        rmsa = wavefield.rmsa(data=self.data, axis=axis, window=window, dim=dim,
                            times=self.times("unix"), distances = self.distances,
                            channels_num=self.channels, vmin=vmin, vmax=vmax,
                            plot_mode=plot_mode)
        if results:
            return rmsa

    def p2p_amp(self, dim="t", results=False, plot_mode="pyqt"):
        """Computes peak-to-peak amplitude of data in time or space.
        See :func:`~fobench.core.fiber.tools.wavefield.peak_to_peak_amp`.
        """
        axis = self.__axis__(dim)
        p2p_amplitude, up_index, down_index = wavefield.peak_to_peak_amp(self.data,
                                             self.sampling_frequency, axis=axis)
        if plot_mode=="pyqt" and dim=="t":
                plot_pyqt.plot_distance(distances=self.distances, channels_num=self.channels,
                            data=p2p_amplitude, y_label="P2P Amplitude", x_label="Channel",
                              title="Peak-to-Peak Amplitude Profile")
        if plot_mode=="pyqt" and dim=="d":
                plot_pyqt.plot_timeseries(timestamps=self.times(time_type="unix"), data=p2p_amplitude,
                                    y_label="Amplitude",
                                    title="Peak-to-Peak Amplitude over time")

        if results:
            return p2p_amplitude, up_index, down_index

    """Plotting methods"""

    def fx_plot(self, norm=False, vmin=None, vmax=None, order=1, nfft=None, figsize=None,
                 show=True, cmap="viridis", results=False, file_name=None,
                 where=None, plot_mode="pyqt", **kwargs):
        """Computes frequency-distance plot.
        See :func:`~fobench.core.tools.wavefield.frequency_content`,
        :func:`~fobench.core.plotting.plotting_mpl.mpl_fx_plot` and
        :func:`~fobench.core.plotting.plotting_pyqt.plot_2d_distance`.
        """
        axis = self.__axis__("t")

        fx, freqs =  wavefield.frequency_content(data=self.data, fs=self.sampling_frequency,
                                           order=order, nfft=nfft, norm=norm, axis=axis)
        p95 = np.percentile(fx, 95)
        if vmin is None: vmin = 0
        if vmax is None: vmax = p95
        if plot_mode == "pyqt":
            plot_pyqt.plot_2d_distance(distances=self.distances, channels_num=np.array(self.channels),
                              y_ticks=freqs, data=fx if axis else fx.T,
                              cmap=cmap, vmin=vmin, vmax=vmax, y_label="Frequency [Hz]",
                              title="Frequency content", cbar_label=self.units)
        elif plot_mode == "mpl":
            plot.mpl_fx_plot(spec_matrix=np.rot90(fx) if axis else fx[::-1], freqs=freqs, x=self.channels,
                     units_y="Energy", figsize=figsize, title=str(self.start_time.date),
                     cmap=cmap, file_name=file_name, vmin=vmin, vmax=vmax, **kwargs)

        if results:
            return fx, freqs

    def spectrum(self, channel, plot_mode="pyqt", norm=False, pre_processing=True,
                 order=1, pad=0, nfft=None, mode="spectrum", figsize=None,
                 nperseg=None, file_name=None, legend=True, results=False, **kwargs):
        """Computes and plots spectrum of channel(s), mode can be ``'spectrum'``
        or ``'psd'``. See :func:`~fobench.core.tools.signals.signal_spectrum`,
        :func:`~fobench.core.plotting.plotting_pyqt.plot_spectral`
        and :func:`~fobench.core.plotting.plotting_mpl.simple_spectrum`.
        """

        axis = self.__axis__("t")
        if isinstance(channel, np.ndarray):
            channel = sorted(channel)
        if isinstance(channel, tuple):
            channel = list(range(min(channel), max(channel) + 1))
        elif isinstance(channel, list):
            channel = sorted(channel)
        else:
            channel = [channel]
        ch_idx = np.array([self.channels.tolist().index(ch) for ch in channel])
        o_signal = np.take(self.data, indices=ch_idx, axis=self.__axis__("d"))

        f, spec = signals.signal_spectrum(o_signal=o_signal, fs=self.sampling_frequency, mode=mode,
                norm=norm, order=order, nfft=nfft, pre_processing=pre_processing, pad=pad, nperseg=nperseg,
                axis=axis)
        spec = spec[:,0] if spec.ndim > 1 else spec  # reduced dimensionality of single spectrum

        if plot_mode=="pyqt":
            units = self.units if mode == "spectrum" else f"{self.units.split(' ')[-1]}²/Hz"
            plot_pyqt.plot_spectral(frequencies=f, amplitudes=spec, y_label =f"{units.title()}",
                                title= f"{mode.title()}" if mode=="spectrum" else f"{mode.upper()}", labels=channel)
        elif plot_mode=="mpl":
            units = self.units if mode == "spectrum" else f"{self.units.split(' ')[-1]}$^{{2}}$/Hz"
            plot.simple_spectrum(spectra=np.array([spec]), freqs=f, channels=[channel], y_units=units, legend=legend, figsize=figsize,
                        title=str(self.start_time.date), file_name=file_name, **kwargs)

        if results:
            return f, spec

    def channel_plot(self, channel, max_value=None, figsize=None, file_name=None,
                    plot_mode="pyqt", **kwargs):
        """Generates simple plot of channel data.
        See :func:`~fobench.core.plotting.plotting_pyqt.plot_timeseries`.
        and :func:`~fobench.core.plotting.plotting_mpl.simple_plot`.
        """
        if isinstance(channel, np.ndarray):
            channel = sorted(channel)
        if isinstance(channel, tuple):
            channel = list(range(min(channel), max(channel) + 1))
        elif isinstance(channel, list):
            channel = sorted(channel)
        else:
            channel = [channel]
        ch_idx = np.array([self.channels.tolist().index(ch) for ch in channel])
        selected = np.take(self.data, indices=ch_idx, axis=self.__axis__("d"))
        if plot_mode=="pyqt":
            t = self.times(time_type="unix")
            plot_pyqt.plot_timeseries(data=selected, timestamps=t, y_label=self.units,
                             dt=self.dt, title="Channel Plot", labels=channel)
        elif plot_mode=="mpl":
            t = self.times("matplotlib")
            plot.simple_plot(data=selected.T, t=t, channel=channel, units_y=self.units,
                    max_value=max_value, figsize=figsize, title=str(self.start_time.date),
                    file_name=file_name, **kwargs)

    def plot(self, vmin=None, vmax=None, figsize=None, show=True, cmap="seismic",
          file_name=None, where=None, add_data=None, plot_mode="pyqt", **kwargs):
        """Generates plot of data. See :func:`~fobench.core.plotting.plotting_pyqt.plot_2d_timeseries`
        and :func:`~fobench.core.plotting.plotting_mpl.gen_DAS_plot`
        """
        if plot_mode == "pyqt":
            t = self.times(time_type="unix")
            p95 = np.percentile(self.data, 95)
            vmin = -p95 if vmin is None else vmin
            vmax = p95 if vmax is None else vmax
            plot_pyqt.plot_2d_timeseries(timestamps=t, y_ticks=np.array(self.channels),
                        data=self.data.T if self.__axis__("t") else self.data,
                        y_label="Channel", dt=self.dt, title="Data Plot", vmin=vmin,
                        vmax=vmax, cbar_label=self.units, distances=self.distances)
        elif plot_mode == "mpl":
            t = self.times(time_type="matplotlib")
            plot.gen_DAS_plot(data=self.data.T if self.__axis__("t") else self.data, t=t,
                        channels=self.channels,     units_y=self.units, figsize=figsize,
                        title=str(self.start_time.date), cmap=cmap, file_name=file_name,
                        vmin=vmin, vmax=vmax, add_data=add_data, **kwargs)

    def channel_spectrogram(self, channel, norm=False, trace=False, figsize=None,
                        cmap="viridis", file_name=None,     freq_lim=None,  results=False,
                        plot_mode="pyqt", vmin=None, vmax=None, **kwargs):
        """Computes and plots spectrogram for a ``"channel"``.
        See :func:`~fobench.core.tools.signals.signal_spectrogram`,
        :func:`~fobench.core.plotting.plotting_pyqt.plot_2d_timeseries` and
        :func:`~fobench.core.plotting.plotting_mpl.simple_spectrogram`.
        """
        axis = self.__axis__("t")
        channel = int(channel)
        index = self.channels.tolist().index(channel)
        data = self.data[:, index]
        f, t, Sxx = signals.signal_spectrogram(data=data, sampling_frequency=self.sampling_frequency,
                                         axis=axis, norm=norm)
        if plot_mode == "pyqt":
            t = self.times(time_type="unix")
            if vmin is None: vmin = 0
            if vmax is None: vmax = np.percentile(Sxx, 95)
            plot_pyqt.plot_2d_timeseries(timestamps=t, y_ticks=f, dt=self.dt,
                        data=np.rot90(Sxx, k=-1), y_label="Frequency [Hz]",
                        title=f"Spectrogram channel {channel}", cmap="viridis",
                        vmin=vmin, vmax=vmax, cbar_label=self.units)
        elif plot_mode == "mpl":
            t = self.times(time_type="matplotlib")
            plot.simple_spectrogram(data=Sxx, freq=f, t=t, units_y=self.units,
                        trace=data if trace == True else None, figsize=figsize, cmap=cmap,
                        title=str(self.start_time.date)+"  "+"Ch:"+str(channel),
                        file_name=file_name, freq_lim=freq_lim, **kwargs)

        if results:
            return Sxx, f, t

    def record_section(self, channels, plot_mode="pyqt"):
        """Plots record section of multiple channels. If channels is tuple the range
        between lower and upper limit will be plotted, if list only channels in list
        will be plotted.
        """
        if isinstance(channels, np.ndarray):
            channels = channels.tolist()
        if isinstance(channels, tuple):
            ch0, chf = sorted(map(int, channels))
            channel_list = self.channels.tolist()
            ch0, chf = channel_list.index(ch0), channel_list.index(chf)
            ch_idx = slice(ch0, chf + 1)
        elif isinstance(channels, list):
            channels = list(map(int, channels))
            channel_list = self.channels.tolist()
            ch_idx = sorted(channel_list.index(ch) for ch in channels)
        else:
            raise TypeError(f"Invalid type for channels: {type(channels).__name__}. Expected tuple, list, or np.ndarray.")
        das_data = self.data[:, ch_idx]
        das_channels = np.array(self.channels)[ch_idx]
        if plot_mode=="pyqt":
            plot_pyqt.plot_record_section(timestamps=self.times("unix"), data=das_data,
                                 title="Record Section", numbers=das_channels, dt=self.dt,
                                 y_label="Channel")
        elif plot_mode=="mpl":
            plot.plot_record_section(signals=das_data, t=self.times("matplotlib"),
                            channels=das_channels, date=str(self.start_time.date))

    def acf_profile(self, max_lag, plot_mode="pyqt", deconvolve=False,
                    window_size=None, results=False, vmin=None, vmax=None, **imshow_kwargs):
        """Computes autocorrelation profile.
        See :func:`~fobench.core.tools.wavefield.autocorrelation_profile`.
        """
        axis = self.__axis__("t")
        max_shift = int(max_lag*self.sampling_frequency)
        if max_shift >= self.num_points:
            raise ValueError("Selected max_shift is too large")
        acf = wavefield.autocorrelation_profile(self.data, max_shift, axis, plot_mode,
                                                deconvolve, self.total_channels,
                                                self.distances, self.channels,
                                                self.sampling_frequency,
                                                window_size=window_size, vmin=vmin,
                                                vmax=vmax, **imshow_kwargs)

        if results:
            return acf

    def spatial_coherence(self, max_lag, results=False, plot_mode="pyqt", vmin=None,
                       vmax=None):
        """Computes sptial coherence matrix.
        See :func:`~fobench.core.tools.wavefield.spatial_coherence_matrix`
        """
        data_input = np.moveaxis(self.data, (self.__axis__("d"), self.__axis__("t")), (0, 1))
        coh = wavefield.spatial_coherence_matrix(data=data_input, max_lag=max_lag,
                                           distances=self.distances,
                                           fs=self.sampling_frequency,
                                           channel_nums=self.channels,
                                           plot_mode=plot_mode, results=results,
                                           vmin=vmin, vmax=vmax)
        if results:
            return coh

    def view(self):
        """Launches the Fobench Data Viewer.
        """
        print(f"{'-'*65}\nStarting Fobench Data Viewer")
        app = QtWidgets.QApplication.instance()
        if app is None:
            app = QtWidgets.QApplication(sys.argv)
        self._viewer = Viewer(self)
        self._viewer.show()
        pg.exec()
        print(f"{'-'*65}")