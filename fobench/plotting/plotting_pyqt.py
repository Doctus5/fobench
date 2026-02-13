'''
contains all functionality related to plotting using PyQtGraph, i.e.
whenever plot_mode is set to 'pyqt'
'''
import sys
import datetime
import numpy as np
from pathlib import Path
from PyQt5 import QtCore
import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets, QtGui

'''
-----------------------------------------------------------------
Line Plot Functions
-----------------------------------------------------------------
'''

def plot_timeseries(timestamps: np.ndarray, data: np.ndarray, dt: float,
                    y_label: str = '', title: str = '') -> None:
    '''
    generate generic time series plot using PyQtGraph, ideal for channel plots

    Parameters
    ----------
    timestamps : np.ndarray
        array containing Unix timestamps of data.
    data : np.ndarray
        array containing data to plot.
    dt : float
        sampling period of data
    y_label : str, optional
        y-axis label. The default is ''.
    title : str, optional
        title of plot. The default is ''.

    Returns
    -------
    None
        -
    '''
    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication(sys.argv)

    win = pg.GraphicsLayoutWidget(show=True)
    win.setWindowTitle(f'Fobench: {title}')
    win.setWindowIcon(QtGui.QIcon(
        str(Path(__file__).resolve().parent / 'logo.png')))
    win.setBackground('w')  # white bg
    win.resize(1200, 500)

    plot = win.addPlot(title=title)
    plot.setTitle(title, size='20pt', color='k')
    plot.setLabel('left', y_label, **{'color': 'k', 'font-size': '14pt'})
    plot.plot(timestamps, data, pen=pg.mkPen('k', width=1))

    # y axis
    y_axis = plot.getAxis('left')
    y_axis.setPen(pg.mkPen('k', width=2))
    y_axis.setStyle(tickFont=pg.Qt.QtGui.QFont('Arial', 14))
    y_axis.setTextPen(pg.mkPen('k'))

    # x-axis
    plot.setAxisItems({'bottom': pg.DateAxisItem()})
    x_axis = plot.getAxis('bottom')
    x_axis.setStyle(tickFont=pg.Qt.QtGui.QFont('Arial', 14))
    x_axis.setPen(pg.mkPen('k', width=2))
    x_axis.setTextPen(pg.mkPen('k'))
    date = datetime.datetime.fromtimestamp(timestamps[0]).strftime('%d.%m.%Y')
    x_axis.setLabel(date, **{'color': 'k', 'font-size': '14pt'})
    plot.setXRange(min(timestamps), max(timestamps), padding=0)
    plot.getViewBox().setLimits(xMin=min(timestamps), xMax=max(timestamps),
                                yMin=min(data), yMax=max(data))

    label = pg.LabelItem(justify='left')
    win.addItem(label, row=2, col=0)

    def mouse_moved(evt, dt=dt):
        pos = evt[0]
        if plot.sceneBoundingRect().contains(pos):
            mouse_point = plot.vb.mapSceneToView(pos)
            x_val = mouse_point.x()
            x_val = round(x_val / dt) * dt
            x_datetime = datetime.datetime.utcfromtimestamp(x_val).strftime('%Y-%m-%d %H:%M:%S.%f')
            label.setText(f'Time: {x_datetime}', color='k')
    proxy = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60, slot=mouse_moved)

    pg.exec()
    QtWidgets.QApplication.processEvents()


def plot_record_section(timestamps: np.ndarray, data: np.ndarray, dt: float,
                        numbers: np.ndarray, y_label: str = '', title: str = '') -> None:
    '''
    extended version of timeseries plot for multi-channel data

    Parameters
    ----------
    timestamps : np.ndarray
        array containing Unix timestamps of data.
    data : np.ndarray
        array containing data to plot.
    dt : float
        sampling period of data
    numbers : np.ndarray
        array containing channel numbers for plotting
    y_label : str, optional
        y-axis label. The default is ''.
    title : str, optional
        title of plot. The default is ''.

    Returns
    -------
    None
        -
    '''

    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication(sys.argv)

    data = data[:, ::-1]
    numbers = numbers[::-1]

    win = pg.GraphicsLayoutWidget(show=True)
    win.setWindowTitle(f'Fobench: {title}')
    win.setWindowIcon(QtGui.QIcon(
        str(Path(__file__).resolve().parent / 'logo.png')))
    win.setBackground('w')  # white bg
    win.resize(1200, 500)

    plot = win.addPlot(title=title)
    plot.setTitle(title, size='20pt', color='k')
    plot.setLabel('left', y_label, **{'color': 'k', 'font-size': '14pt'})

    offset = data[:].max()
    ticks = []
    for i in range(data.shape[1]):
        y = -i * offset
        plot.plot(timestamps, data[:, i] + y, pen=pg.mkPen('k', width=1))
        ticks.append((y, str(numbers[i])))

    # y axis
    y_axis = plot.getAxis('left')
    y_axis.setTicks([ticks])
    y_axis.setPen(pg.mkPen('k', width=2))
    y_axis.setStyle(tickFont=pg.Qt.QtGui.QFont('Arial', 14))
    y_axis.setTextPen(pg.mkPen('k'))
    y_axis.enableAutoSIPrefix(False)
    y_axis.setStyle(showValues=True)

    # x-axis
    plot.setAxisItems({'bottom': pg.DateAxisItem(utcOffset=1)})
    x_axis = plot.getAxis('bottom')
    x_axis.setStyle(tickFont=pg.Qt.QtGui.QFont('Arial', 14))
    x_axis.setPen(pg.mkPen('k', width=2))
    x_axis.setTextPen(pg.mkPen('k'))
    date = datetime.datetime.fromtimestamp(timestamps[0]).strftime('%d.%m.%Y')
    x_axis.setLabel(date, **{'color': 'k', 'font-size': '14pt'})
    plot.setXRange(min(timestamps), max(timestamps), padding=0)
    plot.getViewBox().setLimits(xMin=min(timestamps), xMax=max(timestamps),
                                yMax=max(data[:, 0]),
                                yMin=min(data[:, i] - i*data[:].max()))

    label = pg.LabelItem(justify='left')
    win.addItem(label, row=2, col=0)

    def mouse_moved(evt, dt=dt):
        pos = evt[0]
        if plot.sceneBoundingRect().contains(pos):
            mouse_point = plot.vb.mapSceneToView(pos)
            x_val = mouse_point.x()
            x_val = round(x_val / dt) * dt
            x_datetime = datetime.datetime.utcfromtimestamp(x_val).strftime('%Y-%m-%d %H:%M:%S.%f')
            label.setText(f'Time: {x_datetime}', color='k')
    proxy = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60, slot=mouse_moved)

    pg.exec()


def plot_distance(distances: np.ndarray, channels_num: np.ndarray, data: np.ndarray,
                  y_label: str = '', x_label: str = 'Channel',
                  title: str = '') -> None:
    '''
    generate generic distances series plot

    Parameters
    ----------
    distances : np.ndarray
        array containing optical distances values
    channels_num : np.ndarray
        array containing channel numbers.
    data : np.ndarray
        array containing data to plot.
    y_label : str, optional
        y-axis label. The default is ''.
    x_label : str, optional
        x-axis label. The default is 'Optical Distance [m]'.
    title : str, optional
        title of plot. The default is ''.

    Returns
    -------
    None
        -
    '''

    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication(sys.argv)

    win = pg.GraphicsLayoutWidget(show=True)
    win.setWindowTitle(f'Fobench: {title}')
    win.setWindowIcon(QtGui.QIcon(
        str(Path(__file__).resolve().parent / 'logo.png')))
    win.setBackground('w')  # white bg
    win.resize(1200, 500)

    dx = channels_num[1] - channels_num[0]

    plot = win.addPlot(title=title)
    plot.setTitle(title, size='20pt', color='k')
    plot.setLabel('left', y_label, **{'color': 'k', 'font-size': '14pt'})
    plot.plot(channels_num, data, pen=pg.mkPen('k', width=1))

    # y axis
    y_axis = plot.getAxis('left')
    y_axis.setPen(pg.mkPen('k', width=2))
    y_axis.setStyle(tickFont=pg.Qt.QtGui.QFont('Arial', 14))
    y_axis.setTextPen(pg.mkPen('k'))

    # x-axis
    x_axis = plot.getAxis('bottom')
    x_axis.setStyle(tickFont=pg.Qt.QtGui.QFont('Arial', 14))
    x_axis.setPen(pg.mkPen('k', width=2))
    x_axis.setTextPen(pg.mkPen('k'))
    x_axis.setLabel(x_label, **{'color': 'k', 'font-size': '14pt'})
    plot.setXRange(min(channels_num), max(channels_num), padding=0)
    plot.getViewBox().setLimits(xMin=min(channels_num), xMax=max(channels_num),
                                yMin=min(data), yMax=max(data))

    # horizontal container
    container = QtWidgets.QWidget()
    h_layout = QtWidgets.QHBoxLayout(container)
    h_layout.setContentsMargins(5, 5, 5, 5)

    # label f cursor tracking
    text_label = QtWidgets.QLabel()
    text_label.setStyleSheet('color: black; font-size: 10pt;')
    h_layout.addWidget(text_label)

    # button for axis switching
    button = QtWidgets.QPushButton('Distance')
    button.setToolTip('Switch x Axis between Channel Number and Optical Distance')
    button.setFixedWidth(70)
    h_layout.addWidget(button)

    # cursor tracking in data units
    def mouse_moved(evt):
        pos = evt[0]
        if plot.sceneBoundingRect().contains(pos):
            mouse_point = plot.vb.mapSceneToView(pos)
            x_val = mouse_point.x()
            x_val = round(x_val / dx) * dx
            text_label.setText(f'{x_axis.labelText}: {x_val:.2f}')
    proxy = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60, slot=mouse_moved)

    def switch_axis():
        nonlocal dx
        if button.text() == 'Distance':
            button.setText('Channel')
            plot.setLabel('bottom', 'Optical Distance [m]', **{'color': 'k', 'font-size': '14pt'})
            x_vals = distances
            dx = distances[1]-distances[0]
        else:
            button.setText('Distance')
            plot.setLabel('bottom', 'Channel', **{'color': 'k', 'font-size': '14pt'})
            x_vals = channels_num
            dx = channels_num[1] - channels_num[0]

        plot.setXRange(min(x_vals), max(x_vals), padding=0)
        plot.getViewBox().setLimits(xMin=min(x_vals), xMax=max(x_vals))
        plot.clear()
        plot.plot(x_vals, data, pen=pg.mkPen('k', width=1))
        plot.enableAutoRange()

    button.clicked.connect(switch_axis)

    proxy_container = QtWidgets.QGraphicsProxyWidget()
    proxy_container.setWidget(container)
    win.addItem(proxy_container, row=2, col=0)

    pg.exec()

def plot_spectral(frequencies: np.ndarray, amplitudes: np.ndarray,
                  y_label: str = 'Amplitude', x_label: str = 'Frequency [Hz]',
                  title: str = '') -> None:
    '''
    generate generic amplitude over frequency plot

    Parameters
    ----------
    frequencies : np.ndarray
        array containing frequency values
    amplitudes : np.ndarray
        array containing amplitudes to plot.
    y_label : str, optional
        y-axis label. The default is 'Amplitude'.
    x_label : str, optional
        x-axis label. The default is 'Frequency [Hz]'.
    title : str, optional
        title of plot. The default is ''.

    Returns
    -------
    None
        -
    '''

    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication(sys.argv)

    win = pg.GraphicsLayoutWidget(show=True)
    win.setWindowTitle(f'Fobench: {title}')
    win.setWindowIcon(QtGui.QIcon(
        str(Path(__file__).resolve().parent / 'logo.png')))
    win.setBackground('w')  # white bg
    win.resize(1200, 500)

    dx = frequencies[1] - frequencies[0]

    plot = win.addPlot(title=title)
    plot.setTitle(title, size='20pt', color='k')
    plot.setLabel('left', y_label, **{'color': 'k', 'font-size': '14pt'})
    plot.plot(frequencies, amplitudes, pen=pg.mkPen('k', width=1))

    # y axis
    y_axis = plot.getAxis('left')
    y_axis.setPen(pg.mkPen('k', width=2))
    y_axis.setStyle(tickFont=pg.Qt.QtGui.QFont('Arial', 14))
    y_axis.setTextPen(pg.mkPen('k'))

    # x-axis
    x_axis = plot.getAxis('bottom')
    x_axis.setStyle(tickFont=pg.Qt.QtGui.QFont('Arial', 14))
    x_axis.setPen(pg.mkPen('k', width=2))
    x_axis.setTextPen(pg.mkPen('k'))
    x_axis.setLabel(x_label, **{'color': 'k', 'font-size': '14pt'})
    plot.setXRange(min(frequencies), max(frequencies), padding=0)
    plot.getViewBox().setLimits(xMin=min(frequencies), xMax=max(frequencies),
                                yMin=min(amplitudes), yMax=max(amplitudes))

    # horizontal container
    container = QtWidgets.QWidget()
    h_layout = QtWidgets.QHBoxLayout(container)
    h_layout.setContentsMargins(5, 5, 5, 5)

    # label f cursor tracking
    text_label = QtWidgets.QLabel()
    text_label.setStyleSheet('color: black; font-size: 10pt;')
    h_layout.addWidget(text_label)

    # cursor tracking in data units
    def mouse_moved(evt):
        pos = evt[0]
        if plot.sceneBoundingRect().contains(pos):
            mouse_point = plot.vb.mapSceneToView(pos)
            x_val = mouse_point.x()
            x_val = round(x_val / dx) * dx
            text_label.setText(f'{x_axis.labelText}: {x_val:.2f}')
    proxy = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60, slot=mouse_moved)

    proxy_container = QtWidgets.QGraphicsProxyWidget()
    proxy_container.setWidget(container)
    win.addItem(proxy_container, row=2, col=0)

    pg.exec()


'''
-----------------------------------------------------------------
Matrix Plot Functions
-----------------------------------------------------------------
'''

def plot_2d_timeseries(timestamps: np.ndarray, data: np.ndarray, y_ticks: list,
                       dt: float, max_value: float = None, y_label: str = '',
                       title: str = '', cmap: str = 'seismic',
                       cbar_label: str = '', distances: np.ndarray = None) -> None:
    '''
    generate generic matrix plot where x-axis represents time, e.g.
    waterfall visualisation of data or spectrograms,PSDs, RMSA ...

    Parameters
    ----------
    timestamps : np.ndarray
        array containing Unix timestamps of data.
    data : np.ndarray
        array containing data to plot.
    y_ticks : list
        y-axis tick labels.
    dt : float
        sampling period of data
    max_value : float, optional
        value that sets limits of colorbar. The default is None.
    y_label : str, optional
        y-axis label. The default is ''.
    title : str, optional
        plot title. The default is ''.
    cmap : str, optional
        colormap to use. The default is 'seismic'.
    cbar_label : str, optional
        label of colorbar. The default is ''.
    distances : np.ndarray, optional
        optional distances to use for y-axis

    Returns
    -------
    None
        -.
    '''

    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication(sys.argv)

    win = pg.GraphicsLayoutWidget(show=True)
    win.setWindowTitle(f'Fobench {title}')
    win.setWindowIcon(QtGui.QIcon(
        str(Path(__file__).resolve().parent / 'logo.png')))
    win.setBackground('w')
    win.resize(1200, 800)

    dy = y_ticks[1] - y_ticks[0]

    plot = win.addPlot(title=title)
    plot.setTitle(title, size='20pt', color='k')
    plot.setLabel('left', y_label, **{'color': 'k', 'font-size': '14pt'})

    plot.setCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))
    plot.setMouseEnabled(x=True, y=True)
    plot.getViewBox().setMouseMode(pg.ViewBox.RectMode)  # One-button

    x_min, x_max = timestamps[0], timestamps[-1]
    y_min, y_max = y_ticks[0], y_ticks[-1]

    image = pg.ImageItem()
    plot.addItem(image)
    image.setImage(data)
    image.setRect(x_min, y_min, x_max - x_min, y_max - y_min)

    # y axis
    y_axis = plot.getAxis('left')
    # y_axis.setTicks([[(y, str(y)) for y in y_ticks]])
    y_axis.setPen(pg.mkPen('k', width=2))
    y_axis.setStyle(tickFont=pg.Qt.QtGui.QFont('Arial', 14))
    y_axis.setTextPen(pg.mkPen('k'))

    # x axis
    plot.setAxisItems({'bottom': pg.DateAxisItem(utcOffset=1)})
    x_axis = plot.getAxis('bottom')
    x_axis.setStyle(tickFont=pg.Qt.QtGui.QFont('Arial', 14))
    x_axis.setTextPen(pg.mkPen('k'))
    x_axis.setPen(pg.mkPen('k', width=2))
    # date = datetime.datetime.fromtimestamp(timestamps[0]).strftime('%d.%m.%Y')
    # x_axis.setLabel(date, **{'color': 'k', 'font-size': '14pt'})

    plot.setXRange(min(timestamps), max(timestamps), padding=0)
    plot.setYRange(min(y_ticks), max(y_ticks), padding=0)
    plot.getViewBox().setLimits(xMin=min(timestamps), xMax=max(timestamps),
                                yMin=min(y_ticks), yMax=max(y_ticks))

    # horizontal container
    container = QtWidgets.QWidget()
    h_layout = QtWidgets.QHBoxLayout(container)
    h_layout.setContentsMargins(5, 5, 5, 5)

    # label f cursor tracking
    text_label = QtWidgets.QLabel()
    text_label.setStyleSheet('color: black; font-size: 10pt;')
    h_layout.addWidget(text_label)

    # button for axis switching if y axis is distance
    if distances:
        button = QtWidgets.QPushButton('Distance')
        button.setToolTip('Switch x Axis between Channel Number and Optical Distance')
        button.setFixedWidth(70)
        h_layout.addWidget(button)

        def switch_axis():
            nonlocal dy
            if button.text() == 'Distance':
                button.setText('Channel')
                plot.setLabel('left', 'Optical Distance [m]', **{'color': 'k', 'font-size': '14pt'})
                y_vals = distances
                dy = distances[1] - distances[0]
            else:
                button.setText('Distance')
                plot.setLabel('left', 'Channel', **{'color': 'k', 'font-size': '14pt'})
                y_vals = y_ticks
                dy = y_ticks[1] - y_ticks[0]
            y_min, y_max = y_vals[0], y_vals[-1]
            image.setRect(x_min, y_min, x_max - x_min, y_max - y_min)
            plot.setYRange(min(y_vals), max(y_vals), padding=0)
            plot.getViewBox().setLimits(yMin=min(y_vals), yMax=max(y_vals))
            plot.enableAutoRange()

        button.clicked.connect(switch_axis)

    def mouse_moved(evt):
        pos = evt[0]
        if plot.sceneBoundingRect().contains(pos):
            mouse_point = plot.vb.mapSceneToView(pos)
            x_val = mouse_point.x()
            x_val = round(x_val / dt) * dt
            y_val = mouse_point.y()
            y_val = round(y_val / dy) * dy

            x_datetime = datetime.datetime.utcfromtimestamp(
                x_val).strftime('%Y-%m-%d %H:%M:%S.%f')
            text_label.setText(f'Time: {x_datetime} | {y_axis.labelText}: {y_val:.2f}')
    proxy = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60, slot=mouse_moved)

    values = (-max_value, max_value) if max_value is not None else None
    data_range = np.nanmax(data) - np.nanmin(data)

    cmap = pg.colormap.get(cmap, source='matplotlib')
    bar = pg.ColorBarItem(colorMap=cmap, values=values,
                          label=cbar_label, interactive=True, rounding=0.001*data_range)
    bar.setImageItem(image, insert_in=plot)

    proxy_container = QtWidgets.QGraphicsProxyWidget()
    proxy_container.setWidget(container)
    win.addItem(proxy_container, row=2, col=0)

    pg.exec()

def plot_2d_distance(distances: np.ndarray, channels_num: np.ndarray,
                     data: np.ndarray, y_ticks: list,
                     max_value: float = None, y_label: str = '',
                     x_label: str = 'Channel',
                     title: str = '', cmap: str = 'seismic',
                     cbar_label: str = '', invert_y=False) -> None:
    '''
    generate generic matrix plot where x-axis represents distance

    Parameters
    ----------
    distances : np.ndarray
        array containing optical distances of channels.
    channels_num : np.ndarray
        array containing channel numbers
    data : np.ndarray
        array containing data to plot.
    y_ticks : list
        y-axis tick labels.
    max_value : float, optional
        value that sets limits of colorbar. The default is None.
    y_label : str, optional
        y-axis label. The default is ''.
    x_label : str, optional
        x-axis label. The default is 'Channel'.
    title : str, optional
        title of plot. The default is ''.
    cmap : str, optional
        colormap to use. The default is 'seismic'.
    cbar_label : str, optional
        label of colorbar. The default is ''.
    invert_y : TYPE, optional
        invert y-axis. The default is False.

    Returns
    -------
    None
        -.
    '''

    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication(sys.argv)

    win = pg.GraphicsLayoutWidget(show=True)
    win.setWindowTitle('Fobench')
    win.setWindowIcon(QtGui.QIcon(
        str(Path(__file__).resolve().parent / 'logo.png')))
    win.setBackground('w')
    win.resize(1200, 800)

    dx = channels_num[1] - channels_num[0]
    dy = y_ticks[1] - y_ticks[0]

    plot = win.addPlot(title=title)
    plot.setTitle(title, size='20pt', color='k')
    plot.setLabel('left', y_label, **{'color': 'k', 'font-size': '14pt'})
    plot.setLabel('bottom', x_label, **{'color': 'k', 'font-size': '14pt'})

    plot.setCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))
    plot.setMouseEnabled(x=True, y=True)
    plot.getViewBox().setMouseMode(pg.ViewBox.RectMode)  # One-button
    if invert_y:
        plot.getViewBox().invertY(True)

    x_min, x_max = channels_num[0], channels_num[-1]
    y_min, y_max = y_ticks[0], y_ticks[-1]

    image = pg.ImageItem()
    plot.addItem(image)
    image.setImage(data)
    image.setRect(x_min, y_min, x_max - x_min, y_max - y_min)

    # y axis
    y_axis = plot.getAxis('left')
    # y_axis.setTicks([[(y, str(y)) for y in y_ticks]])
    y_axis.setPen(pg.mkPen('k', width=2))
    y_axis.setStyle(tickFont=pg.Qt.QtGui.QFont('Arial', 14))
    y_axis.setTextPen(pg.mkPen('k'))

    # x axis
    x_axis = plot.getAxis('bottom')
    x_axis.setStyle(tickFont=pg.Qt.QtGui.QFont('Arial', 14))
    x_axis.setTextPen(pg.mkPen('k'))
    x_axis.setPen(pg.mkPen('k', width=2))

    plot.setXRange(min(channels_num), max(channels_num), padding=0)
    plot.setYRange(min(y_ticks), max(y_ticks), padding=0)
    plot.getViewBox().setLimits(xMin=min(channels_num), xMax=max(channels_num),
                                yMin=min(y_ticks), yMax=max(y_ticks))

    # horizontal container
    container = QtWidgets.QWidget()
    h_layout = QtWidgets.QHBoxLayout(container)
    h_layout.setContentsMargins(5, 5, 5, 5)

    # label f cursor tracking
    text_label = QtWidgets.QLabel()
    text_label.setStyleSheet('color: black; font-size: 10pt;')
    h_layout.addWidget(text_label)

    # button for axis switching
    button = QtWidgets.QPushButton('Distance')
    button.setToolTip('Switch x Axis between Channel Number and Optical Distance')
    button.setFixedWidth(70)
    h_layout.addWidget(button)

    # cursor tracking in data units
    def mouse_moved(evt):
        pos = evt[0]
        if plot.sceneBoundingRect().contains(pos):
            mouse_point = plot.vb.mapSceneToView(pos)
            x_val = mouse_point.x()
            x_val = round(x_val / dx) * dx
            y_val = mouse_point.y()
            y_val = round(y_val / dy) * dy
            text_label.setText(f'{x_axis.labelText}: {x_val:.2f} | {y_label}: {y_val:.3f}')

    proxy = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60, slot=mouse_moved)

    values = (-max_value, max_value) if max_value is not None else None
    data_range = np.nanmax(data) - np.nanmin(data)

    cmap = pg.colormap.get(cmap, source='matplotlib')
    bar = pg.ColorBarItem(colorMap=cmap, values=values, label=cbar_label,
                          interactive=True, rounding=0.001*data_range)
    bar.setImageItem(image, insert_in=plot)

    def switch_axis():
        nonlocal dx
        if button.text() == 'Distance':
            button.setText('Channel')
            plot.setLabel(
                'bottom', 'Optical Distance [m]', **{'color': 'k', 'font-size': '14pt'})
            x_vals = distances
            dx = distances[1] - distances[0]
        else:
            button.setText('Distance')
            plot.setLabel('bottom', 'Channel', **{'color': 'k', 'font-size': '14pt'})
            x_vals = channels_num
            dx = channels_num[1] - channels_num[0]
        x_min, x_max = x_vals[0], x_vals[-1]
        image.setRect(x_min, y_min, x_max - x_min, y_max - y_min)
        plot.setXRange(min(x_vals), max(x_vals), padding=0)
        plot.getViewBox().setLimits(xMin=min(x_vals), xMax=max(x_vals))
        plot.enableAutoRange()

    button.clicked.connect(switch_axis)

    # add container
    proxy_container = QtWidgets.QGraphicsProxyWidget()
    proxy_container.setWidget(container)
    win.addItem(proxy_container, row=2, col=0)

    pg.exec()