'''
contains the Explorer class
'''

import datetime
import numpy as np
import pyqtgraph as pg
from pathlib import Path
from pyqtgraph.Qt import QtCore, QtWidgets, QtGui
from pyqtgraph.dockarea.Dock import Dock
from pyqtgraph.dockarea.DockArea import DockArea


class Explorer(QtWidgets.QMainWindow):
    def __init__(self, Fiber):
        super().__init__()

        # set background color to white
        pg.setConfigOption('background', 'w')
        pg.setConfigOption('foreground', 'k')

        # set up Fiber class information
        self.Fiber = Fiber
        self.times = self.Fiber.times(time_type='unix')
        self.ch0, self.chf = self.Fiber.channels[0], self.Fiber.channels[-1]
        self.selected_channel = self.ch0
        self.selected_distance = self.Fiber.distances[0]

        # set up window
        self.setWindowTitle('Fobench Data Explorer')
        self.setGeometry(100, 100, 1200, 800)
        self.setWindowIcon(QtGui.QIcon(str(Path(__file__).resolve().parent / 'logo.png')))

        self.area = DockArea()
        self.setCentralWidget(self.area)

        # set up docks
        # main dock for data view
        self.dock_1 = Dock('Data View', size=(1200,600), closable=False)
        self.area.addDock(self.dock_1)
        # channel view dock
        self.dock_2 = Dock(f'Channel {self.selected_channel}', size=(1200,200), closable=False)
        self.area.addDock(self.dock_2, 'bottom', self.dock_1)
        # controls dock
        self.dock_3 = Dock('Controls', size=(150,200), closable=False)
        self.area.addDock(self.dock_3, 'right', self.dock_2)

        # boundaries for main data plot
        x_min, x_max = self.times[0], self.times[-1]
        y_min, y_max = self.Fiber.distances[0], self.Fiber.distances[-1]

        # main data plot
        self.matrix_plot_widget = pg.GraphicsLayoutWidget()
        self.matrix_plot = self.matrix_plot_widget.addPlot()
        self.matrix_image = pg.ImageItem(image=self.Fiber.data)
        self.matrix_plot.addItem(self.matrix_image)
        self.matrix_plot.setAspectLocked(False)
        self.matrix_plot.setAxisItems({'bottom': pg.DateAxisItem(utcOffset=1)})
        self.matrix_plot.setLabel('left', 'Optical Distance [m]', **{'color': 'k', 'font-size': '14pt'})
        self.matrix_image.setRect(x_min, y_min, x_max - x_min, y_max - y_min)
        self.matrix_plot.getViewBox().setLimits(xMin=x_min, xMax=x_max,
                                                yMin=y_min, yMax=y_max)
        self.dock_1.addWidget(self.matrix_plot_widget)


        # line tool channel selector
        self.y_picker = pg.InfiniteLine(pos=self.Fiber.distances[0], angle=0, movable=True,
                markers=[('^', 0, 16), ('v', 1, 16)], pen=pg.mkPen('k', width=1, style=pg.QtCore.Qt.DashLine),
                hoverPen=pg.mkPen('k', width=1, style=pg.QtCore.Qt.SolidLine))
        self.y_picker.setBounds([y_min, y_max])
        self.matrix_plot.addItem(self.y_picker)
        self.y_picker.sigDragged.connect(self.on_picker_moved)

        # cursor tracking in data units
        self.matrix_label = pg.LabelItem(justify='left')
        self.matrix_plot_widget.addItem(self.matrix_label, row=1, col=0)

        def mouse_moved(evt):
            pos = evt[0]
            mouse_point = self.matrix_plot.getViewBox().mapSceneToView(pos)
            x_val = mouse_point.x()
            y_val = mouse_point.y()
            x_datetime = datetime.datetime.utcfromtimestamp(x_val).strftime('%Y-%m-%d %H:%M:%S.%f')
            self.matrix_label.setText(f'Time: {x_datetime} | Optical Distance [m]: {y_val:.1f}', color='k')

        self.matrix_mouse_proxy = pg.SignalProxy(self.matrix_plot.scene().sigMouseMoved,
            rateLimit=60, slot=mouse_moved)

        # histogram to main plot
        self.hist = pg.HistogramLUTItem(image=self.matrix_image)
        self.hist.setToolTip('Right click to change colormap')
        self.hist.gradient.setColorMap(pg.colormap.get('seismic', source='matplotlib'))
        self.matrix_plot_widget.addItem(self.hist)

        # channel plot
        self.line_plot_widget = pg.GraphicsLayoutWidget()
        self.line_plot = self.line_plot_widget.addPlot()
        self.line_plot.setAxisItems({'bottom': pg.DateAxisItem(utcOffset=1)})
        self.line_curve = self.line_plot.plot(pen='k')
        self.dock_2.addWidget(self.line_plot_widget)

        # channel selection spinbox
        self.spin_box = QtWidgets.QSpinBox()
        self.spin_box.setToolTip('Selects Channel to plot')
        self.spin_box.setMinimum(self.ch0)
        self.spin_box.setMaximum(self.chf)
        self.spin_box.setPrefix('Channel: ')
        self.spin_box.setValue(self.selected_channel)
        self.spin_box.valueChanged.connect(self.on_spin_box_changed)
        self.dock_3.addWidget(self.spin_box, row=0, col=0)
        self.distance_indicator = QtWidgets.QLabel(f'Distance: {self.selected_distance}')
        self.dock_3.addWidget(self.distance_indicator, row=1, col=0)

        # colorbar max percentile spinbox
        self.cbar_spinbox = QtWidgets.QDoubleSpinBox()
        self.cbar_spinbox.setRange(0.0, 100.0)
        self.cbar_spinbox.setDecimals(0)
        self.cbar_spinbox.setSingleStep(1)
        self.cbar_spinbox.setValue(90)  # default to 90th percentile
        self.cbar_spinbox.setPrefix('Colorbar max: ')
        self.cbar_spinbox.setSuffix(' %')
        self.cbar_spinbox.setToolTip('Sets data percentile to be represented by colorbar (0–100)')
        self.cbar_spinbox.valueChanged.connect(self.update_colorbar_levels)
        self.dock_3.addWidget(self.cbar_spinbox, row=2, col=0)
        self.update_colorbar_levels()

        # dropdown channel analysis menu
        self.dropdown_ch = QtWidgets.QComboBox()
        self.dropdown_ch.addItems(['Channel Analysis', 'Spectrogram', 'PSD', 'Spectrum'])
        self.dropdown_ch.setToolTip('Basic Signal Analysis Methods')
        self.dropdown_ch.currentIndexChanged.connect(self.on_dropdown_ch_changed)
        self.update_plots()
        self.dock_3.addWidget(self.dropdown_ch, row=3, col=0)

        # dropdown data analysis menu
        self.dropdown_data = QtWidgets.QComboBox()
        self.dropdown_data.addItems(['Data Analysis', 'RMSA', 'P2PA', 'f-x Plot'])
        self.dropdown_data.setToolTip('Basic Data Analysis Methods')
        self.dropdown_data.currentIndexChanged.connect(self.on_dropdown_data_changed)
        self.update_plots()
        self.dock_3.addWidget(self.dropdown_data, row=4, col=0)

    # detect channel change with moving horizontal line
    def on_picker_moved(self):
        idx = np.abs(np.asarray(self.Fiber.distances) - self.y_picker.getYPos()).argmin()
        self.selected_channel = self.Fiber.channels[idx]
        self.update_plots()

    # detect channel selection with spin box
    def on_spin_box_changed(self, value):
        self.selected_channel = int(value)
        self.update_plots()

    # detect selection of method in dropdown channel menu
    def on_dropdown_ch_changed(self, index):
        if self.dropdown_ch.currentText() == 'Spectrogram':
            self.plot_spectrogram()
        if self.dropdown_ch.currentText() == 'PSD':
            self.plot_psd()
        if self.dropdown_ch.currentText() == 'Spectrum':
            self.plot_spectrum()
        self.dropdown_ch.setCurrentIndex(0)  # Reset dropdown menus

    # detect selection of method in dropdown data menu
    def on_dropdown_data_changed(self, index):
        if self.dropdown_data.currentText() == 'RMSA':
            self.plot_rmsa()
        if self.dropdown_data.currentText() == 'P2PA':
            self.plot_p2pa()
        if self.dropdown_data.currentText() == 'f-x Plot':
            self.plot_fx()
        self.dropdown_data.setCurrentIndex(0)

    # update all plots and widgets
    def update_plots(self):
        ind = self.selected_channel-self.ch0
        self.selected_distance = self.Fiber.distances[ind]
        self.spin_box.setValue(self.selected_channel)
        self.distance_indicator.setText(f'Distance: {self.selected_distance}')
        self.dock_2.setTitle(f'Channel {self.selected_channel}')
        self.y_picker.setValue(self.selected_distance)
        row_data = self.Fiber.data[:, ind]
        self.line_curve.setData(y=row_data, x=self.times)
        x_min, x_max = self.times[0], self.times[-1]
        y_min, y_max = min(row_data), max(row_data)
        self.line_plot.setXRange(self.times[0], self.times[-1], padding=0)
        self.line_plot.getViewBox().setLimits(xMin=x_min, xMax=x_max,
                                                yMin=y_min, yMax=y_max)

    def plot_spectrogram(self):
        spectrogram_data, freqs = self.spectrogram(self.selected_channel)
        spectrogram_data = np.flip(spectrogram_data, axis=1)

        spectrogram_plot_widget = pg.GraphicsLayoutWidget()
        spectrogram_plot = spectrogram_plot_widget.addPlot()
        spectrogram_image = pg.ImageItem()
        spectrogram_plot.addItem(spectrogram_image)
        spectrogram_plot.setAspectLocked(False)
        spectrogram_plot.setLabel('left', 'Frequency [Hz]')
        spectrogram_plot.setAxisItems({'bottom': pg.DateAxisItem(utcOffset=1)})
        spectrogram_image.setImage(spectrogram_data, levels=(spectrogram_data.min(), spectrogram_data.max()))
        lut = pg.colormap.get('viridis', source='matplotlib').getLookupTable()
        spectrogram_image.setLookupTable(lut)

        x_min, x_max = self.times[0], self.times[-1]
        y_min, y_max = min(freqs), max(freqs)

        spectrogram_image.setRect(x_min, y_min, x_max - x_min, y_max - y_min)
        spectrogram_plot.getViewBox().setLimits(xMin=x_min, xMax=x_max,
                                                yMin=y_min, yMax=y_max)

        self.dock_spec = Dock(f'Spectrogram Ch {self.selected_channel}', size=(1200,200), closable=True)
        self.area.addDock(self.dock_spec, 'above', self.dock_2)
        self.dock_spec.addWidget(spectrogram_plot_widget)

    def plot_fx(self):
        fx, freqs = self.fx()

        fx_plot_widget = pg.GraphicsLayoutWidget()
        fx_plot = fx_plot_widget.addPlot()
        fx_image = pg.ImageItem()
        fx_plot.addItem(fx_image)
        fx_plot.setAspectLocked(False)
        fx_plot.setLabel('left', 'Frequency [Hz]')
        fx_plot.setLabel('bottom', 'Optical Distance [m]')
        fx_image.setImage(fx.T, levels=(fx.min(), fx.max()))
        lut = pg.colormap.get('viridis', source='matplotlib').getLookupTable()
        fx_image.setLookupTable(lut)

        x_min, x_max = self.Fiber.distances[0], self.Fiber.distances[-1]
        y_min, y_max = min(freqs), max(freqs)

        fx_image.setRect(x_min, y_min, x_max - x_min, y_max - y_min)
        fx_plot.getViewBox().setLimits(xMin=x_min, xMax=x_max,
                                                yMin=y_min, yMax=y_max)
        data_range = np.nanmax(fx) - np.nanmin(fx)
        cmap = pg.colormap.get('viridis', source='matplotlib')
        bar = pg.ColorBarItem(colorMap=cmap, interactive=True, rounding=0.001*data_range)
        bar.setImageItem(fx_image, insert_in=fx_plot)

        self.dock_fx = Dock(f'Frequency-Distance', size=(1200,600), closable=True)
        self.area.addDock(self.dock_fx, 'above', self.dock_1)
        self.dock_fx.addWidget(fx_plot_widget)

    def plot_rmsa(self):
        rmsa_data = self.Fiber.rmsa(results=True, plot_mode=None)[1][0,:]

        rmsa_plot_widget = pg.GraphicsLayoutWidget()
        rmsa_plot = rmsa_plot_widget.addPlot()
        rmsa_plot.setLabel('bottom', 'Optical Distance [m]', **{'color': 'k', 'font-size': '10pt'})
        rmsa_plot.setLabel('left', 'RMS Amplitude', **{'color': 'k', 'font-size': '10pt'})

        rmsa_curve = rmsa_plot.plot(self.Fiber.distances, rmsa_data, pen='k')
        rmsa_plot.setXRange(self.Fiber.distances[0], self.Fiber.distances[-1], padding=0)
        rmsa_plot.getViewBox().setLimits(xMin=self.Fiber.distances[0], xMax=self.Fiber.distances[-1],
                                    yMin=min(rmsa_data), yMax=max(rmsa_data))

        self.dock_rmsa = Dock('RMS Amplitude', size=(1200,200), closable=True)
        self.area.addDock(self.dock_rmsa, 'above', self.dock_2)
        self.dock_rmsa.addWidget(rmsa_plot_widget)


    def plot_p2pa(self):
        p2pa_plot_widget = pg.GraphicsLayoutWidget()
        p2pa_plot = p2pa_plot_widget.addPlot()
        p2pa_plot.setLabel('bottom', 'Optical Distance [m]', **{'color': 'k', 'font-size': '10pt'})
        p2pa_plot.setLabel('left', 'P2P Amplitude', **{'color': 'k', 'font-size': '10pt'})

        p2pa_data = self.Fiber.p2p_amp(plot_mode=None)[0]

        p2pa_curve = p2pa_plot.plot(self.Fiber.distances, p2pa_data, pen='k')
        p2pa_plot.setXRange(self.Fiber.distances[0], self.Fiber.distances[-1], padding=0)
        p2pa_plot.getViewBox().setLimits(xMin=self.Fiber.distances[0], xMax=self.Fiber.distances[-1],
                                    yMin=min(p2pa_data), yMax=max(p2pa_data))

        self.dock_p2p = Dock('P2P Amplitude', size=(1200,200), closable=True)
        self.area.addDock(self.dock_p2p, 'above', self.dock_2)
        self.dock_p2p.addWidget(p2pa_plot_widget)


    def plot_psd(self):
        psd_plot_widget = pg.GraphicsLayoutWidget()
        psd_plot = psd_plot_widget.addPlot()
        psd_plot.setLabel('bottom', 'Frequency [Hz]', **{'color': 'k', 'font-size': '10pt'})
        psd_plot.setLabel('left', 'PSD Amplitude', **{'color': 'k', 'font-size': '10pt'})

        freqs, amps = self.Fiber.spectrum(self.selected_channel, mode='psd',
                                          plot_mode=None, results=True)

        psd_curve = psd_plot.plot(freqs, amps, pen='k')
        psd_plot.setXRange(freqs[0], freqs[-1], padding=0)
        psd_plot.getViewBox().setLimits(xMin=freqs[0], xMax=freqs[-1],
                                    yMin=min(amps), yMax=max(amps))

        self.dock_psd = Dock(f'PSD Ch {self.selected_channel}', size=(1200,200), closable=True)
        self.area.addDock(self.dock_psd, 'above', self.dock_2)
        self.dock_psd.addWidget(psd_plot_widget)

    def plot_spectrum(self):
        spec_plot_widget = pg.GraphicsLayoutWidget()
        spec_plot = spec_plot_widget.addPlot()
        spec_plot.setLabel('bottom', 'Frequency [Hz]', **{'color': 'k', 'font-size': '10pt'})
        spec_plot.setLabel('left', 'FFT Amplitude', **{'color': 'k', 'font-size': '10pt'})

        freqs, amps = self.Fiber.spectrum(self.selected_channel, mode='spectrum',
                                          plot_mode=None, results=True)

        spec_curve = spec_plot.plot(freqs, amps, pen='k')
        spec_plot.setXRange(freqs[0], freqs[-1], padding=0)
        spec_plot.getViewBox().setLimits(xMin=freqs[0], xMax=freqs[-1],
                                    yMin=min(amps), yMax=max(amps))

        self.dock_fft = Dock(f'FFT Ch {self.selected_channel}', size=(1200,200), closable=True)
        self.area.addDock(self.dock_fft, 'above', self.dock_2)
        self.dock_fft.addWidget(spec_plot_widget)

    def spectrogram(self, channel):
        Sxx, f, t = self.Fiber.channel_spectrogram(channel, results=True, make_plot=False)
        return Sxx.T, f

    def fx(self):
        return self.Fiber.fx_plot(plot_mode=None, results=True)

    def update_colorbar_levels(self):
        percentile = self.cbar_spinbox.value()
        data = self.Fiber.data
        abs_percentile = np.percentile(np.abs(data), percentile)
        self.hist.setLevels(-abs_percentile, abs_percentile)