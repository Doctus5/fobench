'''
Contains the Explorer class
'''

import datetime
import numpy as np
import pyqtgraph as pg
from pathlib import Path
from pyqtgraph.Qt import QtWidgets, QtGui
from pyqtgraph.dockarea.Dock import Dock
from pyqtgraph.dockarea.DockArea import DockArea


class Explorer(QtWidgets.QMainWindow):
    '''
    A GUI for exploring FOS data
    '''
    def __init__(self, Fiber):
        super().__init__()

        # set background color to white
        pg.setConfigOptions(background='w', foreground='k')

        # set up Fiber class information
        self.Fiber = Fiber
        self.times = Fiber.times(time_type='unix')
        self.ch0, self.chf = Fiber.channels[0], Fiber.channels[-1]
        self.selected_channel = self.ch0
        self.selected_distance = Fiber.distances[0]

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
        y_min, y_max = Fiber.distances[0], Fiber.distances[-1]

        # main data plot
        self.matrix_plot_widget = pg.GraphicsLayoutWidget()
        self.matrix_plot = self.matrix_plot_widget.addPlot()
        self.matrix_image = pg.ImageItem(image=Fiber.data)
        self.matrix_image.setRect(x_min, y_min, x_max - x_min, y_max - y_min)
        self.matrix_plot.addItem(self.matrix_image)
        self.matrix_plot.setAspectLocked(False)
        self.matrix_plot.setAxisItems({'bottom': pg.DateAxisItem(utcOffset=1)})
        self.matrix_plot.setLabel('left', 'Optical Distance [m]', color='k', font_size='14pt')
        self.matrix_plot.getViewBox().setLimits(xMin=x_min, xMax=x_max,
                                                yMin=y_min, yMax=y_max)
        self.dock_1.addWidget(self.matrix_plot_widget)

        # line tool channel selector
        self.y_picker = pg.InfiniteLine(pos=Fiber.distances[0], angle=0, movable=True,
                markers=[('^', 0, 16), ('v', 1, 16)], pen=pg.mkPen('k', width=1, style=pg.QtCore.Qt.DashLine),
                hoverPen=pg.mkPen('k', width=1, style=pg.QtCore.Qt.SolidLine))
        self.y_picker.setBounds([y_min, y_max])
        self.y_picker.sigDragged.connect(self.on_picker_moved)
        self.matrix_plot.addItem(self.y_picker)

        # cursor tracking in data units
        self.matrix_label = pg.LabelItem(justify='left')
        self.matrix_plot_widget.addItem(self.matrix_label, row=1, col=0)

        def mouse_moved(evt):
            mouse_point = self.matrix_plot.getViewBox().mapSceneToView(evt[0])
            x_datetime = datetime.datetime.utcfromtimestamp(mouse_point.x()).strftime('%Y-%m-%d %H:%M:%S.%f')
            self.matrix_label.setText(f'Time: {x_datetime} | Optical Distance [m]: {mouse_point.y():.1f}', color='k')

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
        self.spin_box.setRange(self.ch0, self.chf)
        self.spin_box.setPrefix('Channel: ')
        self.spin_box.setValue(self.selected_channel)
        self.spin_box.valueChanged.connect(self.on_spin_box_changed)
        self.dock_3.addWidget(self.spin_box, row=0, col=0)
        self.distance_indicator = QtWidgets.QLabel(f'Distance: {self.selected_distance}')
        self.dock_3.addWidget(self.distance_indicator, row=1, col=0)

        # colorbar max percentile spinbox
        self.cbar_spinbox = QtWidgets.QSpinBox()
        self.cbar_spinbox.setRange(0, 100)
        self.cbar_spinbox.setValue(95)  # default to 90th percentile
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
        self.dock_3.addWidget(self.dropdown_ch, row=3, col=0)

        # dropdown data analysis menu
        self.dropdown_data = QtWidgets.QComboBox()
        self.dropdown_data.addItems(['Data Analysis', 'RMSA', 'P2PA', 'f-x Plot'])
        self.dropdown_data.setToolTip('Basic Data Analysis Methods')
        self.dropdown_data.currentIndexChanged.connect(self.on_dropdown_data_changed)
        self.dock_3.addWidget(self.dropdown_data, row=4, col=0)

        self.update_plots()

    def on_picker_moved(self):
        '''Detects selection of channel selection via horizontal line'''
        idx = np.abs(np.asarray(self.Fiber.distances) - self.y_picker.getYPos()).argmin()
        self.selected_channel = self.Fiber.channels[idx]
        self.update_plots()

    # detect
    def on_spin_box_changed(self, value):
        '''Detects selection of channel selection via spin box'''
        self.selected_channel = int(value)
        self.update_plots()

    def on_dropdown_ch_changed(self, index):
        '''Detects selection of method in dropdown channel menu'''
        actions = {'Spectrogram': self.plot_spectrogram,
                   'PSD': self.plot_psd,
                   'Spectrum': self.plot_spectrum,}
        action = actions.get(self.dropdown_ch.currentText())
        if action: action()
        self.dropdown_ch.setCurrentIndex(0)

    #
    def on_dropdown_data_changed(self, index):
        '''Detects selection of method in dropdown data menu'''
        actions = {'RMSA': self.plot_rmsa,
                   'P2PA': self.plot_p2pa,
                   'f-x Plot': self.plot_fx,}
        action = actions.get(self.dropdown_data.currentText())
        if action: action()
        self.dropdown_data.setCurrentIndex(0)

    def update_plots(self):
        '''Updates all plots'''
        ind = self.selected_channel-self.ch0
        self.selected_distance = self.Fiber.distances[ind]
        self.spin_box.setValue(self.selected_channel)
        self.distance_indicator.setText(f'Distance: {self.selected_distance}')
        self.dock_2.setTitle(f'Channel {self.selected_channel}')
        self.y_picker.setValue(self.selected_distance)
        row_data = self.Fiber.data[:, ind]
        self.line_curve.setData(x=self.times, y=row_data)
        self.line_plot.setXRange(self.times[0], self.times[-1], padding=0)
        self.line_plot.getViewBox().setLimits(
            xMin=self.times[0], xMax=self.times[-1],
            yMin=row_data.min(), yMax=row_data.max())

    def plot_spectrogram(self):
        '''Simple spectrogram plot'''
        spectrogram_data, freqs = self.spectrogram(self.selected_channel)
        spectrogram_data = np.flip(spectrogram_data, axis=1)
        x_min, x_max = self.times[0], self.times[-1]
        y_min, y_max = freqs.min(), freqs.max()

        spectrogram_plot_widget = pg.GraphicsLayoutWidget()
        spectrogram_plot = spectrogram_plot_widget.addPlot()
        spectrogram_plot.setAspectLocked(False)
        spectrogram_plot.setLabel('left', 'Frequency [Hz]')
        spectrogram_plot.setAxisItems({'bottom': pg.DateAxisItem(utcOffset=1)})
        spectrogram_image = pg.ImageItem()
        spectrogram_image.setImage(spectrogram_data, levels=(spectrogram_data.min(), spectrogram_data.max()))
        spectrogram_image.setLookupTable(pg.colormap.get('viridis', source='matplotlib').getLookupTable())
        spectrogram_image.setRect(x_min, y_min, x_max - x_min, y_max - y_min)
        spectrogram_plot.addItem(spectrogram_image)
        spectrogram_plot.getViewBox().setLimits(xMin=x_min, xMax=x_max, yMin=y_min, yMax=y_max)
        self.dock_spec = Dock(f'Spectrogram Ch {self.selected_channel}', size=(1200, 200), closable=True)
        self.area.addDock(self.dock_spec, 'above', self.dock_2)
        self.dock_spec.addWidget(spectrogram_plot_widget)

    def plot_fx(self):
        '''Frequency-Distance plot'''
        fx, freqs = self.Fiber.fx_plot(plot_mode=None, results=True)
        x_min, x_max = self.Fiber.distances[0], self.Fiber.distances[-1]
        y_min, y_max = freqs.min(), freqs.max()

        fx_plot_widget = pg.GraphicsLayoutWidget()
        fx_plot = fx_plot_widget.addPlot()
        fx_plot.setAspectLocked(False)
        fx_plot.setLabel('left', 'Frequency [Hz]')
        fx_plot.setLabel('bottom', 'Optical Distance [m]')
        cmap = pg.colormap.get('viridis', source='matplotlib')
        fx_image = pg.ImageItem()
        fx_image.setImage(fx.T, levels=(fx.min(), fx.max()))
        fx_image.setLookupTable(cmap.getLookupTable())
        fx_image.setRect(x_min, y_min, x_max - x_min, y_max - y_min)
        fx_plot.addItem(fx_image)
        bar = pg.ColorBarItem(colorMap=cmap, interactive=True, rounding=0.001 * (np.nanmax(fx) - np.nanmin(fx)))
        bar.setImageItem(fx_image, insert_in=fx_plot)
        fx_plot.getViewBox().setLimits(xMin=x_min, xMax=x_max, yMin=y_min, yMax=y_max)
        self.dock_fx = Dock('Frequency-Distance', size=(1200, 600), closable=True)
        self.area.addDock(self.dock_fx, 'above', self.dock_1)
        self.dock_fx.addWidget(fx_plot_widget)

    def plot_rmsa(self):
        '''RMS amplitude plot'''
        rmsa_data = self.Fiber.rmsa(results=True, plot_mode=None)[1][0, :]
        x_min, x_max = self.Fiber.distances[0], self.Fiber.distances[-1]
        rmsa_plot_widget = pg.GraphicsLayoutWidget()
        rmsa_plot = rmsa_plot_widget.addPlot()
        rmsa_plot.setLabel('bottom', 'Optical Distance [m]', color='k', font_size='10pt')
        rmsa_plot.setLabel('left', 'RMS Amplitude', color='k', font_size='10pt')
        rmsa_plot.plot(self.Fiber.distances, rmsa_data, pen='k')
        rmsa_plot.setXRange(x_min, x_max, padding=0)
        rmsa_plot.getViewBox().setLimits(xMin=x_min, xMax=x_max, yMin=rmsa_data.min(), yMax=rmsa_data.max())
        self.dock_rmsa = Dock('RMS Amplitude', size=(1200, 200), closable=True)
        self.area.addDock(self.dock_rmsa, 'above', self.dock_2)
        self.dock_rmsa.addWidget(rmsa_plot_widget)

    def plot_p2pa(self):
        '''Peak-to-Peak amplitude plot'''
        p2pa_plot_widget = pg.GraphicsLayoutWidget()
        p2pa_plot = p2pa_plot_widget.addPlot()
        p2pa_plot.setLabel('bottom', 'Optical Distance [m]', color='k', font_size='10pt')
        p2pa_plot.setLabel('left', 'P2P Amplitude', color='k', font_size='10pt')
        p2pa_data = self.Fiber.p2p_amp(plot_mode=None)[0]
        p2pa_plot.plot(self.Fiber.distances, p2pa_data, pen='k')
        p2pa_plot.setXRange(self.Fiber.distances[0], self.Fiber.distances[-1], padding=0)
        p2pa_plot.getViewBox().setLimits(xMin=self.Fiber.distances[0], xMax=self.Fiber.distances[-1],
                                    yMin=min(p2pa_data), yMax=max(p2pa_data))
        self.dock_p2p = Dock('P2P Amplitude', size=(1200,200), closable=True)
        self.area.addDock(self.dock_p2p, 'above', self.dock_2)
        self.dock_p2p.addWidget(p2pa_plot_widget)


    def plot_psd(self):
        '''Power spectral density plot'''
        psd_plot_widget = pg.GraphicsLayoutWidget()
        psd_plot = psd_plot_widget.addPlot()
        psd_plot.setLabel('bottom', 'Frequency [Hz]', color='k', font_size='10pt')
        psd_plot.setLabel('left', 'PSD Amplitude', color='k', font_size='10pt')
        freqs, amps = self.Fiber.spectrum(self.selected_channel, mode='psd',
                                          plot_mode=None, results=True)
        psd_plot.plot(freqs, amps, pen='k')
        psd_plot.setXRange(freqs[0], freqs[-1], padding=0)
        psd_plot.getViewBox().setLimits(xMin=freqs[0], xMax=freqs[-1],
                                    yMin=min(amps), yMax=max(amps))
        self.dock_psd = Dock(f'PSD Ch {self.selected_channel}', size=(1200,200), closable=True)
        self.area.addDock(self.dock_psd, 'above', self.dock_2)
        self.dock_psd.addWidget(psd_plot_widget)

    def plot_spectrum(self):
        '''Amplitude spectrum plot'''
        spec_plot_widget = pg.GraphicsLayoutWidget()
        spec_plot = spec_plot_widget.addPlot()
        spec_plot.setLabel('bottom', 'Frequency [Hz]', color='k', font_size='10pt')
        spec_plot.setLabel('left', 'FFT Amplitude', color='k', font_size='10pt')
        freqs, amps = self.Fiber.spectrum(self.selected_channel, mode='spectrum',
                                          plot_mode=None, results=True)
        spec_plot.plot(freqs, amps, pen='k')
        spec_plot.setXRange(freqs[0], freqs[-1], padding=0)
        spec_plot.getViewBox().setLimits(xMin=freqs[0], xMax=freqs[-1],
                                    yMin=min(amps), yMax=max(amps))
        self.dock_fft = Dock(f'FFT Ch {self.selected_channel}', size=(1200,200), closable=True)
        self.area.addDock(self.dock_fft, 'above', self.dock_2)
        self.dock_fft.addWidget(spec_plot_widget)

    def spectrogram(self, channel):
        '''Calls Fiber.spectrogram and returns the result'''
        Sxx, f, _ = self.Fiber.channel_spectrogram(channel, results=True, make_plot=False)
        return Sxx.T, f

    def update_colorbar_levels(self):
        '''Updates the colorbar levels after change by user'''
        percentile = self.cbar_spinbox.value()
        abs_percentile = np.percentile(np.abs(self.Fiber.data), percentile)
        self.hist.setLevels(-abs_percentile, abs_percentile)