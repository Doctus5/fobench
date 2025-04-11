'''
contains the Explorer class
'''

import sys
sys.path.append('/home/joni/Dokumente/GEO4D/Software/fobench/')

from pathlib import Path
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtWidgets, QtGui

class Explorer(QtWidgets.QMainWindow):
    def __init__(self, Fiber):
        super().__init__()
        
        # Set background color to white
        pg.setConfigOption('background', 'w')
        pg.setConfigOption('foreground', 'k')
        
        # set up Fiber class information
        self.Fiber = Fiber
        self.times = self.Fiber.times(time_type='unix')
        self.selected_channel = 0
        
        x_min, x_max = self.times[0], self.times[-1]
        y_min, y_max = self.Fiber.distances[0], self.Fiber.distances[-1]
        
        # set up window
        self.setWindowTitle('Fobench Data Explorer')
        self.setGeometry(100, 100, 1200, 800)
        self.setWindowIcon(QtGui.QIcon(str(Path(__file__).resolve().parent / 'logo.png')))

        central_widget = QtWidgets.QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QtWidgets.QVBoxLayout(central_widget)

        # main data plot
        self.matrix_plot_widget = pg.GraphicsLayoutWidget()
        main_layout.addWidget(self.matrix_plot_widget, stretch=3)  # matrix plot is 3 parts of the layout
        self.matrix_plot = self.matrix_plot_widget.addPlot()
        self.matrix_image = pg.ImageItem(image=self.Fiber.data)
        self.matrix_plot.addItem(self.matrix_image)
        self.matrix_plot.setAspectLocked(False)
        self.matrix_plot.scene().sigMouseClicked.connect(self.on_matrix_clicked)
        self.matrix_plot.setAxisItems({'bottom': pg.DateAxisItem(utcOffset=1)}) # when is utcoffset necessary??
        self.matrix_plot.setLabel('left', 'Optical Distance [m]', **{'color': 'k', 'font-size': '14pt'})
        self.matrix_image.setRect(x_min, y_min, x_max - x_min, y_max - y_min)
        self.matrix_plot.getViewBox().setLimits(xMin=x_min, xMax=x_max,
                                                yMin=y_min, yMax=y_max)

        # histogram to main plot
        self.hist = pg.HistogramLUTItem(image=self.matrix_image)
        self.hist.gradient.setColorMap(pg.colormap.get('seismic', source='matplotlib'))
        self.matrix_plot_widget.addItem(self.hist)
        level_val = max(abs(self.Fiber.data.min()), self.Fiber.data.max())

        # TODO: improve colorbar values 
        self.hist.setLevels(0.1 * -level_val, 0.1 * level_val)
        
        # tabs
        self.tab_widget = QtWidgets.QTabWidget()
        main_layout.addWidget(self.tab_widget, stretch=1)  # Tab widget takes 1 part of the layout

        self.line_plot_widget = pg.GraphicsLayoutWidget()
        self.line_plot = self.line_plot_widget.addPlot()
        self.line_plot.setAxisItems({'bottom': pg.DateAxisItem(utcOffset=1)})
        
        self.line_curve = self.line_plot.plot(pen='k')
        self.tab_widget.addTab(self.line_plot_widget, 'Channel')

        bottom_layout = QtWidgets.QHBoxLayout()
        main_layout.addLayout(bottom_layout)

        self.channel_label = QtWidgets.QLabel('Channel:')
        bottom_layout.addWidget(self.channel_label)  # Channel label in the same line as spinbox and slider

        self.spin_box = QtWidgets.QSpinBox()
        self.spin_box.setMinimum(0)
        self.spin_box.setMaximum(self.Fiber.data.shape[1] - 1)
        self.spin_box.setValue(self.selected_channel)
        self.spin_box.valueChanged.connect(self.on_spin_box_changed)
        bottom_layout.addWidget(self.spin_box)  # Spin box to the left of the slider

        self.slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(self.Fiber.data.shape[1] - 1)
        self.slider.setValue(self.selected_channel)
        self.slider.valueChanged.connect(self.on_slider_changed)
        bottom_layout.addWidget(self.slider)  # Slider to the right of the spin box

        self.dropdown = QtWidgets.QComboBox()
        self.dropdown.addItems(['Methods', 'Spectrogram', 'PSD', 'RMSA', 'P2PA'])
        self.dropdown.currentIndexChanged.connect(self.on_dropdown_changed)
        bottom_layout.addWidget(self.dropdown)  # Dropdown menu to the right of the slider

        self.update_plots()

    def on_matrix_clicked(self, event):
        if event.button() == QtCore.Qt.LeftButton:
            pos = self.matrix_image.mapFromScene(event.scenePos())
            x, y = int(pos.x()), int(pos.y())
            if 0 <= y < self.Fiber.data.shape[1]:
                self.selected_channel = y
                self.slider.setValue(self.selected_channel)
                self.spin_box.setValue(self.selected_channel)
                self.update_plots()

    def on_slider_changed(self, value):
        self.selected_channel = value
        self.spin_box.setValue(self.selected_channel)
        self.update_plots()

    def on_spin_box_changed(self, value):
        self.selected_channel = int(value)
        self.slider.setValue(self.selected_channel)
        self.update_plots()

    def on_dropdown_changed(self, index):
        if self.dropdown.currentText() == 'Spectrogram':
            self.plot_spectrogram()
        if self.dropdown.currentText() == 'RMSA':
            self.plot_rmsa()
        if self.dropdown.currentText() == 'P2PA':
            self.plot_p2pa()
        if self.dropdown.currentText() == 'PSD':
            self.plot_psd()
        self.dropdown.setCurrentIndex(0)  # Reset dropdown to 'Methods'

    def update_plots(self):
        row_data = self.Fiber.data[:, self.selected_channel]
        self.line_curve.setData(y=row_data, x=self.times)
        x_min, x_max = self.times[0], self.times[-1]
        y_min, y_max = min(row_data), max(row_data)
        self.line_plot.setXRange(self.times[0], self.times[-1], padding=0)
        self.line_plot.getViewBox().setLimits(xMin=x_min, xMax=x_max,
                                                yMin=y_min, yMax=y_max)

    def plot_spectrogram(self):
        spectrogram_data, freqs = self.spectrogram(self.selected_channel)
        spectrogram_data = np.flip(spectrogram_data, axis=1)
        # Check if the spectrogram tab already exists
        for i in range(self.tab_widget.count()):
            if self.tab_widget.tabText(i) == f'Spectrogram {self.selected_channel}':
                self.tab_widget.setCurrentIndex(i)
                return

        # Create a new tab for the spectrogram plot
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

        # Add a close button to the tab
        close_button = QtWidgets.QToolButton()
        close_button.setText('x')
        close_button.clicked.connect(lambda: self.close_tab(self.tab_widget.indexOf(spectrogram_plot_widget)))

        self.tab_widget.addTab(spectrogram_plot_widget, f'Spectrogram Ch {self.selected_channel}')
        self.tab_widget.tabBar().setTabButton(self.tab_widget.count() - 1, QtWidgets.QTabBar.RightSide, close_button)
        self.tab_widget.setCurrentWidget(spectrogram_plot_widget)

    def plot_rmsa(self):
        # Check if the RMSA tab already exists
        for i in range(self.tab_widget.count()):
            if self.tab_widget.tabText(i) == 'RMSA':
                self.tab_widget.setCurrentIndex(i)
                return

        # Create a new tab for the RMSA plot
        rmsa_plot_widget = pg.GraphicsLayoutWidget()
        rmsa_plot = rmsa_plot_widget.addPlot()
        rmsa_plot.setLabel('bottom', 'Optical Distance [m]', **{'color': 'k', 'font-size': '10pt'})
        rmsa_plot.setLabel('left', 'RMS Amplitude', **{'color': 'k', 'font-size': '10pt'})

        # Get the RMSA data
        rmsa_data = self.Fiber.rmsa()[1][0]

        # Plot the RMSA data
        rmsa_curve = rmsa_plot.plot(self.Fiber.distances, rmsa_data, pen='k')
        rmsa_plot.setXRange(self.Fiber.distances[0], self.Fiber.distances[-1], padding=0)
        rmsa_plot.getViewBox().setLimits(xMin=self.Fiber.distances[0], xMax=self.Fiber.distances[-1],
                                    yMin=min(rmsa_data), yMax=max(rmsa_data))

        # Add a close button to the tab
        close_button = QtWidgets.QToolButton()
        close_button.setText('x')
        close_button.clicked.connect(lambda: self.close_tab(self.tab_widget.indexOf(rmsa_plot_widget)))

        # Add the new tab to the tab widget
        self.tab_widget.addTab(rmsa_plot_widget, 'RMSA')
        self.tab_widget.tabBar().setTabButton(self.tab_widget.count() - 1, QtWidgets.QTabBar.RightSide, close_button)
        self.tab_widget.setCurrentWidget(rmsa_plot_widget)

    def plot_p2pa(self):
        # Check if the tab already exists
        for i in range(self.tab_widget.count()):
            if self.tab_widget.tabText(i) == 'P2PA':
                self.tab_widget.setCurrentIndex(i)
                return

        # Create a new tab for the p2pa plot
        p2pa_plot_widget = pg.GraphicsLayoutWidget()
        p2pa_plot = p2pa_plot_widget.addPlot()

        # Get the p2pa data
        p2pa_data = self.Fiber.pp_amp()

        # Plot the p2pa data
        p2pa_curve = p2pa_plot.plot(self.Fiber.distances, p2pa_data, pen='k')
        p2pa_plot.setLabel('bottom', 'Optical Distance [m]', **{'color': 'k', 'font-size': '10pt'})
        p2pa_plot.setLabel('left', 'P2P Amplitude', **{'color': 'k', 'font-size': '10pt'})
        
        p2pa_plot.setXRange(self.Fiber.distances[0], self.Fiber.distances[-1], padding=0)
        p2pa_plot.getViewBox().setLimits(xMin=self.Fiber.distances[0], xMax=self.Fiber.distances[-1],
                                    yMin=min(p2pa_data), yMax=max(p2pa_data))


        # Add a close button to the tab
        close_button = QtWidgets.QToolButton()
        close_button.setText('X')
        close_button.clicked.connect(lambda: self.close_tab(self.tab_widget.indexOf(p2pa_plot_widget)))

        # Add the new tab to the tab widget
        self.tab_widget.addTab(p2pa_plot_widget, 'P2PA')
        self.tab_widget.tabBar().setTabButton(self.tab_widget.count() - 1, QtWidgets.QTabBar.RightSide, close_button)
        self.tab_widget.setCurrentWidget(p2pa_plot_widget)
        
    def plot_psd(self):
        # Create a new tab for the p2pa plot
        psd_plot_widget = pg.GraphicsLayoutWidget()
        psd_plot = psd_plot_widget.addPlot()

        freqs, amps = self.Fiber.spectrum(self.selected_channel, s_type='psd', results=True)
        
        psd_curve = psd_plot.plot(freqs, amps, pen='k')
        psd_plot.setLabel('bottom', 'Frequency [Hz]', **{'color': 'k', 'font-size': '10pt'})
        psd_plot.setLabel('left', 'PSD Amplitude', **{'color': 'k', 'font-size': '10pt'})
        
        psd_plot.setXRange(freqs[0], freqs[-1], padding=0)
        psd_plot.getViewBox().setLimits(xMin=freqs[0], xMax=freqs[-1],
                                    yMin=min(amps), yMax=max(amps))

        # Add a close button to the tab
        close_button = QtWidgets.QToolButton()
        close_button.setText('X')
        close_button.clicked.connect(lambda: self.close_tab(self.tab_widget.indexOf(psd_plot_widget)))

        # Add the new tab to the tab widget
        self.tab_widget.addTab(psd_plot_widget, f'PSD Ch: {self.selected_channel}')
        self.tab_widget.tabBar().setTabButton(self.tab_widget.count() - 1, QtWidgets.QTabBar.RightSide, close_button)
        self.tab_widget.setCurrentWidget(psd_plot_widget)
            



    def close_tab(self, index):
        self.tab_widget.removeTab(index)

    def spectrogram(self, channel):
        Sxx, f, t = self.Fiber.channel_spectrogram(channel, verbose=True, make_plot=False)
        return Sxx.T, f

