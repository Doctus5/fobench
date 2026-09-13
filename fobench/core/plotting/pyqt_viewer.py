"""Contains the Viewer class."""

import datetime
import functools
from pathlib import Path
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets, QtGui
from pyqtgraph.dockarea import Dock, DockArea

def _busy_cursor(func):
    """Decorator to call methods with pg.BusyCursor."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        with pg.BusyCursor():
            return func(*args, **kwargs)
    return wrapper


class Viewer(QtWidgets.QMainWindow):
    """Interactive GUI for exploring fiber optic sensing data."""
    def __init__(self, Fiber: object) -> None:
        super().__init__()
        pg.setConfigOptions(background="w", foreground="k")

        # set up Fiber class information
        self.Fiber = Fiber
        self.times = Fiber.times(time_type="unix")
        self.ch0, self.chf = Fiber.channels[0], Fiber.channels[-1]
        self.selected_channel = self.ch0
        self.selected_distance = Fiber.distances[0]
        self.data = Fiber.data.T if Fiber.__axis__("t") else Fiber.data

        # set up window
        self.setWindowTitle("Fobench Data Viewer")
        self.setWindowIcon(QtGui.QIcon(str(Path(__file__).resolve().parent/"logo.png")))

        # set up the docks
        self.area = DockArea()
        self.setCentralWidget(self.area)
        self.dock_3 = Dock("Controls", size=(80, 800), closable=False)
        self.area.addDock(self.dock_3, "right")
        self.dock_1 = Dock("Data View", size=(1000, 600), closable=False)
        self.area.addDock(self.dock_1, "left", self.dock_3)
        self.dock_2 = Dock(f"Channel {self.selected_channel}", size=(1000, 200), closable=False)
        self.area.addDock(self.dock_2, "bottom", self.dock_1)

        # boundaries for main data plot
        x_min, x_max = self.times[0], self.times[-1]
        y_min, y_max = Fiber.distances[0], Fiber.distances[-1]

        # state of y-axis
        self._y_state = {"y_vals": np.array(Fiber.distances, dtype=float),
                         "dy": Fiber.distances[1]-Fiber.distances[0]}

        # main data plot
        self.matrix_plot_widget = pg.GraphicsLayoutWidget()
        self.matrix_plot = self.matrix_plot_widget.addPlot()
        self.matrix_image = pg.ImageItem(image=self.data)
        self.matrix_image.setRect(x_min, y_min, x_max-x_min, y_max-y_min)
        self.matrix_plot.addItem(self.matrix_image)
        self.matrix_plot.setAspectLocked(False)
        self.matrix_plot.setAxisItems({"bottom": pg.DateAxisItem(utcOffset=1)})
        self.matrix_plot.setLabel("left", "Optical Distance [m]", color="k", font_size="14pt")
        self.matrix_plot.getViewBox().setLimits(xMin=x_min, xMax=x_max,
                                                yMin=y_min, yMax=y_max)
        self.dock_1.addWidget(self.matrix_plot_widget)

        # line tool channel selector
        self.y_picker = pg.InfiniteLine(pos=Fiber.distances[0], angle=0, movable=True,
                markers=[("^", 0, 16), ("v", 1, 16)], pen=pg.mkPen("k", width=1, style=pg.QtCore.Qt.DashLine),
                hoverPen=pg.mkPen("k", width=1, style=pg.QtCore.Qt.SolidLine))
        self.y_picker.setBounds([y_min, y_max])
        self.y_picker.sigDragged.connect(self.on_picker_moved)
        self.matrix_plot.addItem(self.y_picker)

        # cursor tracking in data units
        self.matrix_label = pg.LabelItem(justify="left")
        self.matrix_plot_widget.addItem(self.matrix_label, row=1, col=0)

        def mouse_moved(evt):
            mouse_point = self.matrix_plot.vb.mapSceneToView(evt[0])
            dt = Fiber.dt
            dy = self._y_state["dy"]
            y_label_text = self.matrix_plot.getAxis("left").labelText
            x_snapped = round(mouse_point.x()/dt)*dt
            x_str = datetime.datetime.utcfromtimestamp(x_snapped).strftime("%Y-%m-%d %H:%M:%S.%f")
            y_snapped = round(mouse_point.y()/dy)*dy
            mouse_point2 = self.matrix_plot.vb.mapFromViewToItem(self.matrix_image, mouse_point)
            col = round(mouse_point2.x())
            row = round(mouse_point2.y())
            if 0 <= col < self.matrix_image.image.shape[0] and 0 <= row < self.matrix_image.image.shape[1]:
                z = self.matrix_image.image[col, row]
                self.matrix_label.setText(
                    f"Time: {x_str} | {y_label_text}: {y_snapped:.1f} | {self.Fiber.units.title()}: {z:e}", color="k")
        self.matrix_mouse_proxy = pg.SignalProxy(self.matrix_plot.scene().sigMouseMoved,
            rateLimit=60, slot=mouse_moved)

        # channel plot
        self.line_plot_widget = pg.GraphicsLayoutWidget()
        self.line_plot = self.line_plot_widget.addPlot()
        self.line_plot.setAxisItems({"bottom": pg.DateAxisItem(utcOffset=1)})
        self.line_curve = self.line_plot.plot(pen="k")
        self.line_plot.setXLink(self.matrix_plot)
        self.line_plot.setLabel("left", f"{self.Fiber.units.title()}")
        self.dock_2.addWidget(self.line_plot_widget)

        # cursor tracking for channel plot
        self.line_label = pg.LabelItem(justify="left")
        self.line_plot_widget.addItem(self.line_label, row=1, col=0)

        def line_mouse_moved(evt):
            mouse_point = self.line_plot.getViewBox().mapSceneToView(evt[0])
            dt = Fiber.dt
            x_snapped = round(mouse_point.x()/dt)*dt
            y = mouse_point.y()
            x_str = datetime.datetime.utcfromtimestamp(x_snapped).strftime("%Y-%m-%d %H:%M:%S.%f")
            self.line_label.setText(f"Time: {x_str} | {self.Fiber.units.title()}: {y:e}", color="k")

        self.line_mouse_proxy = pg.SignalProxy(self.line_plot.scene().sigMouseMoved,
            rateLimit=60, slot=line_mouse_moved)

        # fill the controls dock
        # histogram
        self.hist = pg.HistogramLUTItem(image=self.matrix_image)
        self.hist.setToolTip("Right click to change colormap")
        self.hist.gradient.setColorMap(pg.colormap.get("seismic", source="matplotlib"))
        self.hist.axis.setLabel(Fiber.units.title())
        self.hist_widget = pg.GraphicsLayoutWidget()
        self.hist_widget.addItem(self.hist)
        self.dock_3.addWidget(self.hist_widget, row=0, col=0)

        # channel selection spinbox
        self.spin_box = QtWidgets.QSpinBox()
        self.spin_box.setToolTip("Selects Channel to plot")
        self.spin_box.setRange(self.ch0, self.chf)
        self.spin_box.setPrefix("Channel: ")
        self.spin_box.setValue(self.selected_channel)
        self.spin_box.valueChanged.connect(self.on_spin_box_changed)
        self.dock_3.addWidget(self.spin_box, row=1, col=0)

        # distance indicator
        self.distance_indicator = QtWidgets.QLabel(f"Distance: {self.selected_distance}")
        self.dock_3.addWidget(self.distance_indicator, row=2, col=0)

        # colorbar max percentile spinbox
        self.cbar_spinbox = QtWidgets.QSpinBox()
        self.cbar_spinbox.setRange(0, 100)
        self.cbar_spinbox.setValue(95)
        self.cbar_spinbox.setPrefix("Colorbar max: ")
        self.cbar_spinbox.setSuffix(" %")
        self.cbar_spinbox.setToolTip("Sets data percentile to be represented by colorbar (0–100)")
        self.cbar_spinbox.valueChanged.connect(self.update_colorbar_levels)
        self.dock_3.addWidget(self.cbar_spinbox, row=3, col=0)

        # dropdown channel analysis menu
        self.dropdown_ch = QtWidgets.QComboBox()
        self.dropdown_ch.addItems(["Channel Analysis", "Spectrogram", "PSD", "Spectrum"])
        self.dropdown_ch.setToolTip("Basic Signal Analysis Methods")
        self.dropdown_ch.currentIndexChanged.connect(self.on_dropdown_ch_changed)
        self.dock_3.addWidget(self.dropdown_ch, row=4, col=0)

        # dropdown data analysis menu
        self.dropdown_data = QtWidgets.QComboBox()
        self.dropdown_data.addItems(["Data Analysis", "RMSA", "P2PA", "SNR", "fx-Plot", "ACFs"])
        self.dropdown_data.setToolTip("Basic Data Analysis Methods")
        self.dropdown_data.currentIndexChanged.connect(self.on_dropdown_data_changed)
        self.dock_3.addWidget(self.dropdown_data, row=5, col=0)

        # dropdown measure tools
        self.dropdown_measure = QtWidgets.QComboBox()
        self.dropdown_measure.addItems(["Measure", "Time-Space", "Velocity"])
        self.dropdown_measure.setToolTip("Measurement Tools")
        self.dropdown_measure.currentIndexChanged.connect(self.on_dropdown_measure_changed)
        self.dock_3.addWidget(self.dropdown_measure, row=6, col=0)

        # link-unlink button
        self.link_button = QtWidgets.QPushButton("Unlink X-Axis")
        self.link_button.setToolTip("Toggle X-axis linking between main and channel plot")
        self.link_button.setCheckable(True)
        self.link_button.clicked.connect(self.on_link_button_clicked)
        self.dock_3.addWidget(self.link_button, row=7, col=0)

        # y-axis switch button
        self.y_axis_button = QtWidgets.QPushButton("Channel")
        self.y_axis_button.setToolTip("Switch Y-axis between Optical Distance and Channel Number")
        self.y_axis_button.clicked.connect(self.on_y_axis_button_clicked)
        self.dock_3.addWidget(self.y_axis_button, row=8, col=0)

        # make enough space on the left
        self.matrix_plot.getAxis("left").setWidth(70)
        self.line_plot.getAxis("left").setWidth(70)

        # update everything and force ranges
        self.update_colorbar_levels()
        self.update_plots()
        self.showMaximized()
        self.matrix_plot.setXRange(self.times[0], self.times[-1], padding=0)
        self.line_plot.setXRange(self.times[0], self.times[-1], padding=0)

    def refresh_y_axis(self, y_vals: np.ndarray) -> None:
        """Handles the y-axis and reposition of matrix plot when changing between distance and channels."""
        x_min, x_max = self.times[0], self.times[-1]
        dy = y_vals[1]-y_vals[0]
        y_min, y_max = y_vals[0], y_vals[-1]
        self._y_state["y_vals"] = y_vals
        self._y_state["dy"] = dy
        self.matrix_image.setRect(x_min, y_min-dy/2, x_max-x_min, y_max-y_min+dy)
        self.matrix_plot.getViewBox().setLimits(xMin=x_min, xMax=x_max,
                                                yMin=-np.inf, yMax=np.inf)
        self.matrix_plot.setYRange(y_min-dy/2, y_max+dy/2, padding=0)
        self.matrix_plot.getViewBox().setLimits(xMin=x_min, xMax=x_max,
                                                yMin=y_min-dy/2, yMax=y_max+dy/2)
        self.y_picker.setBounds([y_min, y_max])

    def on_y_axis_button_clicked(self) -> None:
        """Switches the main plot y-axis between optical distance and channel number."""
        if self.y_axis_button.text() == "Channel":
            self.y_axis_button.setText("Distance")
            self.matrix_plot.setLabel("left", "Channel", color="k", font_size="14pt")
            self.refresh_y_axis(np.array(self.Fiber.channels, dtype=float))
            self.y_picker.setValue(float(self.selected_channel))
        else:
            self.y_axis_button.setText("Channel")
            self.matrix_plot.setLabel("left", "Optical Distance [m]", color="k", font_size="14pt")
            self.refresh_y_axis(np.array(self.Fiber.distances, dtype=float))
            self.y_picker.setValue(self.selected_distance)

    def on_picker_moved(self) -> None:
        """Detects selection of channel via horizontal line."""
        if self.y_axis_button.text() == "Distance":
            idx = np.abs(np.array(self.Fiber.channels, dtype=float)-self.y_picker.getYPos()).argmin()
        else:
            idx = np.abs(np.asarray(self.Fiber.distances)-self.y_picker.getYPos()).argmin()
        self.selected_channel = self.Fiber.channels[idx]
        self.update_plots()

    def on_spin_box_changed(self, value: float) -> None:
        """Detects selection of channel selection via spin box."""
        self.selected_channel = int(value)
        self.update_plots()

    def on_dropdown_ch_changed(self, index: int) -> None:
        """Detects selection of method in dropdown channel menu."""
        actions = {"Spectrogram": self.plot_spectrogram,
                   "PSD": lambda: self.plot_spectrum(mode="psd"),
                   "Spectrum": lambda: self.plot_spectrum(mode="spectrum")}
        action = actions.get(self.dropdown_ch.currentText())
        if action: action()
        self.dropdown_ch.setCurrentIndex(0)

    def on_dropdown_data_changed(self, index:int) -> None:
        """Detects selection of method in dropdown data menu."""
        actions = {"RMSA":  lambda: self.plot_data_over_distance(mode="rmsa"),
                   "P2PA":  lambda: self.plot_data_over_distance(mode="p2p_amp"),
                   "SNR":  lambda: self.plot_data_over_distance(mode="snr"),
                   "fx-Plot": self.plot_fx,
                   "ACFs": self.acf_dialog}
        action = actions.get(self.dropdown_data.currentText())
        if action: action()
        self.dropdown_data.setCurrentIndex(0)

    def on_dropdown_measure_changed(self, index:int) -> None:
        """Detects selection of method in dropdown measure menu."""
        actions = {"Time-Space":  self.add_rect_roi,
                   "Velocity": self.add_velocity_line}
        action = actions.get(self.dropdown_measure.currentText())
        self.dropdown_measure.setCurrentIndex(0)
        if action: action()

    def update_plots(self) -> None:
        """Updates all plots."""
        ind = self.selected_channel-self.ch0
        self.selected_distance = self.Fiber.distances[ind]
        self.spin_box.blockSignals(True)
        self.spin_box.setValue(self.selected_channel)
        self.spin_box.blockSignals(False)
        self.distance_indicator.setText(f"Distance: {self.selected_distance}")
        self.dock_2.setTitle(f"Channel {self.selected_channel}")
        if self.y_axis_button.text() == "Distance":
            self.y_picker.setValue(float(self.selected_channel))
        else:
            self.y_picker.setValue(self.selected_distance)
        row_data = self.data[:, ind]
        self.line_curve.setData(x=self.times, y=row_data)
        self.line_plot.getViewBox().setLimits(
            xMin=self.times[0], xMax=self.times[-1],
            yMin=row_data.min(), yMax=row_data.max())

    @_busy_cursor
    def plot_spectrogram(self) -> None:
        """Simple spectrogram plot."""
        spectrogram_data, freqs, _ = self.Fiber.channel_spectrogram(self.selected_channel,
                                                                    results=True, plot_mode=None)
        spectrogram_data = np.flip(spectrogram_data.T, axis=1)
        x_min, x_max = self.times[0], self.times[-1]
        y_min, y_max = freqs.min(), freqs.max()
        cmap = pg.colormap.get("viridis", source="matplotlib")
        spec_max = np.percentile(spectrogram_data, 95)
        spectrogram_plot_widget = pg.GraphicsLayoutWidget()
        spectrogram_plot = spectrogram_plot_widget.addPlot()
        spectrogram_plot.setAspectLocked(False)
        spectrogram_plot.setLabel("left", "Frequency [Hz]")
        spectrogram_plot.setAxisItems({"bottom": pg.DateAxisItem(utcOffset=1)})
        spectrogram_image = pg.ImageItem()
        spectrogram_image.setImage(spectrogram_data, levels=(0, spec_max))
        spectrogram_image.setLookupTable(cmap.getLookupTable())
        spectrogram_image.setRect(x_min, y_min, x_max-x_min, y_max-y_min)
        spectrogram_plot.addItem(spectrogram_image)
        bar = pg.ColorBarItem(colorMap=cmap, values=(0, spec_max),
                              interactive=True, label=self.Fiber.units.title(),
                              rounding=0.001*(spectrogram_data.max()-spectrogram_data.min()))
        bar.setImageItem(spectrogram_image, insert_in=spectrogram_plot)
        spectrogram_plot.getViewBox().setLimits(xMin=x_min, xMax=x_max, yMin=y_min, yMax=y_max)
        self.dock_spec = Dock(f"Spectrogram Ch {self.selected_channel}", size=(1200, 200), closable=True)
        self.area.addDock(self.dock_spec, "above", self.dock_2)
        self.dock_spec.addWidget(spectrogram_plot_widget)

    @_busy_cursor
    def plot_fx(self) -> None:
        """Frequency-Distance plot."""
        fx, freqs = self.Fiber.fx_plot(plot_mode=None, results=True)
        x_min, x_max = self.Fiber.distances[0], self.Fiber.distances[-1]
        y_min, y_max = freqs.min(), freqs.max()

        fx_plot_widget = pg.GraphicsLayoutWidget()
        fx_plot = fx_plot_widget.addPlot()
        fx_plot.setAspectLocked(False)
        fx_plot.setLabel("left", "Frequency [Hz]")
        fx_plot.setLabel("bottom", "Optical Distance [m]")
        cmap = pg.colormap.get("viridis", source="matplotlib")
        fx_image = pg.ImageItem()
        fx_abs_max = np.percentile(fx, 95)
        fx_image.setImage(fx.T)
        fx_image.setLookupTable(cmap.getLookupTable())
        fx_image.setRect(x_min, y_min, x_max-x_min, y_max-y_min)
        fx_plot.addItem(fx_image)
        fx_min, fx_max = float(np.nanmin(fx)), float(np.nanmax(fx))
        bar = pg.ColorBarItem(colorMap=cmap, values=(0, fx_abs_max),
                              interactive=True, rounding=0.001*(fx_max-fx_min),
                              label=self.Fiber.units.title())
        bar.setImageItem(fx_image, insert_in=fx_plot)
        fx_plot.getViewBox().setLimits(xMin=x_min, xMax=x_max, yMin=y_min, yMax=y_max)
        self.dock_fx = Dock("Frequency-Distance", size=(1200, 600), closable=True)
        self.area.addDock(self.dock_fx, "above", self.dock_1)
        self.dock_fx.addWidget(fx_plot_widget)

    def acf_dialog(self) -> None:
        "Plot the Autocorrelation profile"
        @_busy_cursor
        def calc(max_lag_val, deconvolve, window_size):
            return self.Fiber.acf_profile(max_lag=max_lag_value,
                                          deconvolve=deconvolve,
                                          window_size=window_size,
                                          plot_mode=None, results=True)
        dialog = QtWidgets.QDialog(self)
        dialog.setWindowTitle("Autocorrelation Settings")
        dialog.setFixedWidth(400)
        layout = QtWidgets.QFormLayout(dialog)
        layout.setSpacing(10)

        max_lag = QtWidgets.QLineEdit("")
        layout.addRow("Maximum Lag (s):", max_lag)
        deconvolve = QtWidgets.QCheckBox("Deconvolve Source term?")
        win_size = QtWidgets.QLineEdit("Full Cable")
        win_size.setEnabled(False)  # start disabled
        layout.addRow(deconvolve)
        layout.addRow("Window size (number of channels)", win_size)
        deconvolve.toggled.connect(win_size.setEnabled)
        buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addRow(buttons)
        if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            max_lag_value = float(max_lag.text())
            window = None if win_size.text() == "Full Cable" else int(win_size.text())
            acfs = calc(max_lag_value, deconvolve.isChecked(), window)
            x_min, x_max = self.Fiber.distances[0], self.Fiber.distances[-1]
            y_min, y_max = 0, max_lag_value
            acf_plot_widget = pg.GraphicsLayoutWidget()
            acf_plot = acf_plot_widget.addPlot()
            acf_plot.setAspectLocked(False)
            acf_plot.setLabel("left", "Lag/TWT [s]")
            acf_plot.setLabel("bottom", "Optical Distance [m]")
            cmap = pg.colormap.get("seismic", source="matplotlib")
            acf_image = pg.ImageItem()
            acf_image.setImage(acfs.T)
            acf_image.setLookupTable(cmap.getLookupTable())
            acf_image.setRect(x_min, y_min, x_max-x_min, y_max-y_min)
            acf_plot.addItem(acf_image)
            bar = pg.ColorBarItem(colorMap=cmap, values=(-0.8, 0.8),
                                  interactive=True, rounding=0.001,
                                  label="Correlation Coefficient")
            bar.setImageItem(acf_image, insert_in=acf_plot)
            acf_plot.getViewBox().setLimits(xMin=x_min, xMax=x_max, yMin=y_min, yMax=y_max)
            acf_plot.getViewBox().invertY(True)
            self.dock_acf = Dock("Autocorrelation Profile", size=(1200, 600), closable=True)
            self.area.addDock(self.dock_acf, "above", self.dock_1)
            self.dock_acf.addWidget(acf_plot_widget)

    def make_plot_dock(self, title: str, x_label: str, y_label: str,
                       dock_ref: Dock = None) -> tuple[pg.PlotItem, Dock]:
        """Returns dock with plot widget."""
        widget = pg.GraphicsLayoutWidget()
        plot = widget.addPlot()
        plot.setLabel("bottom", x_label, color="k", font_size="10pt")
        plot.setLabel("left", y_label, color="k", font_size="10pt")
        dock = Dock(title, size=(1200, 200), closable=True)
        self.area.addDock(dock, "above", dock_ref or self.dock_2)
        dock.addWidget(widget)
        return plot, dock

    @_busy_cursor
    def plot_data_over_distance(self, mode) -> None:
        """Plots P2PA-, RMSA-, and SNR-profile."""
        if mode == "rmsa":
            data = self.Fiber.rmsa(results=True, plot_mode=None)[0, :]
            plot, self.dock_rmsa = self.make_plot_dock("RMS Amplitude", "Optical Distance [m]",
                                                       "RMS Amplitude")
        elif mode == "p2p_amp":
            data, _, _ = self.Fiber.p2p_amp(plot_mode=None, results=True)
            plot, self.dock_p2p = self.make_plot_dock("P2P Amplitude", "Optical Distance [m]",
                                                       "P2P Amplitude")
        elif mode == "snr":
            data = self.Fiber.SNR(plot_mode=None, results=True)
            plot, self.dock_snr = self.make_plot_dock("SNR", "Optical Distance [m]",
                                                       "SNR")

        x_min, x_max = self.Fiber.distances[0], self.Fiber.distances[-1]
        plot.plot(self.Fiber.distances, data, pen="k")
        plot.setXRange(x_min, x_max, padding=0)
        plot.getViewBox().setLimits(xMin=x_min, xMax=x_max,
                                    yMin=min(data), yMax=max(data))

    def plot_spectrum(self, mode: str = None) -> None:
        """Amplitude spectrum or PSD plot."""
        freqs, amps = self.Fiber.spectrum(self.selected_channel, mode=mode,
                                          plot_mode=None, results=True)
        if mode == "spectrum":
            plot, self.dock_fft = self.make_plot_dock(f"FFT Ch {self.selected_channel}",
                                                       "Frequency [Hz]", "FFT Amplitude")
        elif mode == "psd":
            plot, self.dock_psd = self.make_plot_dock(f"PSD Ch {self.selected_channel}",
                                                       "Frequency [Hz]", "PSD Amplitude")
        plot.plot(freqs, np.squeeze(amps), pen="k")
        plot.setXRange(freqs[0], freqs[-1], padding=0)
        plot.getViewBox().setLimits(xMin=freqs[0], xMax=freqs[-1],
                                    yMin=min(amps), yMax=max(amps))

    def update_colorbar_levels(self) -> None:
        """Updates the colorbar levels after change by user."""
        percentile = self.cbar_spinbox.value()
        abs_percentile = np.percentile(np.abs(self.data), percentile)
        self.hist.setLevels(-abs_percentile, abs_percentile)


    def on_link_button_clicked(self) -> None:
        """Handles (un)linking of x-axis."""
        if self.link_button.isChecked():
            self.line_plot.setXLink(None)
            self.link_button.setText("Link X-Axis")
        else:
            self.line_plot.setXLink(self.matrix_plot)
            self.link_button.setText("Unlink X-Axis")

    def add_rect_roi(self) -> pg.RectROI:
        """Adds a rectangular ROI to the matrix plot to measure time and distance
        and displays its time range on the line plot.
        """

        roi = pg.RectROI([self.times[0], self.Fiber.distances[0]],  # x, y
                         [(self.times[-1] - self.times[0]) * 0.1,       # width
                          (self.Fiber.distances[-1] - self.Fiber.distances[0]) * 0.1],  # height
                         pen=pg.mkPen("k", width=2),
                         handlePen=pg.mkPen("k", width=2),
                         hoverPen=pg.mkPen("k", width=2, style=pg.QtCore.Qt.DashLine,),
                         movable=True, resizable=True, removable=True,
                         maxBounds=pg.QtCore.QRectF(self.times[0],
                                                    self.Fiber.distances[0],
                                                    self.times[-1] - self.times[0],
                                                    self.Fiber.distances[-1] - self.Fiber.distances[0]))
        roi.addScaleHandle([0,0], center=[1,1])
        self.matrix_plot.addItem(roi)

        label = pg.TextItem(color="k", anchor=(0, 0),
                            fill=pg.mkBrush(255, 255, 255, 200))
        self.matrix_plot.addItem(label)

        line_roi = pg.LinearRegionItem(values=[self.times[0],
                self.times[0] + (self.times[-1] - self.times[0]) * 0.1],
                                       orientation="vertical",
                                       brush=pg.mkBrush(0, 0, 0, 40),
                                       pen=pg.mkPen("k", width=1), movable=False)
        self.line_plot.addItem(line_roi)

        def update_roi():
            """Updates ROI label and secondary ROI"""
            pos, size = roi.pos(), roi.size()
            x0 = pos.x()
            x1 = x0 + size.x()
            label.setText(f"Duration: {self.format_time(size.x())}\n"
                          f"Distance: {self.format_distance(size.y())}")
            label.setPos(pos.x(), pos.y() + size.y())
            line_roi.setRegion([x0, x1])

        def remove_roi():
            """Handles removing everything upon context menu click"""
            self.matrix_plot.removeItem(roi)
            self.matrix_plot.removeItem(label)
            self.line_plot.removeItem(line_roi)

        roi.sigRegionChanged.connect(update_roi)
        update_roi()
        roi.sigRemoveRequested.connect(remove_roi)

    def add_velocity_line(self) -> pg.LineSegmentROI:
        """Add a line to the matrix plot to measure velocity."""
        x0, y0 = self.times[0], self.Fiber.distances[0]
        x1 = x0 + (self.times[-1] - self.times[0]) * 0.1
        y1 = y0 + (self.Fiber.distances[-1] - self.Fiber.distances[0]) * 0.1

        line = pg.LineSegmentROI([[x0, y0], [x1, y1]],
                                 pen=pg.mkPen("k", width=2),
                                 handlePen=pg.mkPen("k", width=2),
                                 hoverPen=pg.mkPen("k", width=2, style=pg.QtCore.Qt.DashLine),
                                 movable=True, removable=True)
        self.matrix_plot.addItem(line)
        label = pg.TextItem(color="k", anchor=(0, 1),
                            fill=pg.mkBrush(255, 255, 255, 200))
        self.matrix_plot.addItem(label)

        def update_label():
            """Updates label that displays velocity"""
            handles = line.getSceneHandlePositions()
            p1, p2 = handles[0][1], handles[1][1]
            p1 = self.matrix_plot.vb.mapSceneToView(p1)
            p2 = self.matrix_plot.vb.mapSceneToView(p2)
            dx = p2.x() - p1.x()
            dy = p2.y() - p1.y()
            if abs(dx) < 1e-12:
                velocity = float("inf")
            else:
                velocity = abs(dy/dx)
            label.setText(f"Velocity: {velocity:.5g} m/s\n"
                          f"Δt: {self.format_time(abs(dx))}\n"
                          f"Δd: {self.format_distance(abs(dy))}")

            label.setPos((p1.x() + p2.x()) / 2,
                         (p1.y() + p2.y()) / 2)

        line.sigRegionChanged.connect(update_label)
        update_label()

        def remove_line():
            """Handles removing everything upon context menu click"""
            self.matrix_plot.removeItem(line)
            self.matrix_plot.removeItem(label)
        line.sigRemoveRequested.connect(remove_line)

    def format_time(self, seconds):
        """Formats time string"""
        if seconds < 1: return f"{seconds * 1000:.3g} ms"
        if seconds < 60: return f"{seconds:.3g} s"
        if seconds < 3600:
            minutes = int(seconds // 60)
            secs = seconds % 60
            if secs < 0.01:
                return f"{minutes} min"
            return f"{minutes} min {secs:.3g} s"
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        if minutes == 0:
            return f"{hours} h"
        return f"{hours} h {minutes} min"

    def format_distance(self, meters):
        """Formats distance string"""
        if meters < 1000: return f"{meters:.3g} m"
        return f"{meters/1000:.3g} km"