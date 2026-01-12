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

def plot_timeseries(timestamps: np.ndarray, data: np.ndarray, y_label: str ='',
                    title:str = '') -> None:
    '''
    generate generic time series plot using PyQtGraph, ideal for channel plots

    Parameters
    ----------
    timestamps : np.ndarray
        array containing Unix timestamps of data.
    data : np.ndarray
        array containing data to plot.
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
    win.setWindowIcon(QtGui.QIcon(str(Path(__file__).resolve().parent / 'logo.png')))
    win.setBackground('w') # white bg
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

    pg.exec()

def plot_distance(distances: np.ndarray, data: np.ndarray, y_label: str ='',
                    x_label: str = '', title:str = '') -> None:
    '''
    generate generic distances series plot using PyQtGraph

    Parameters
    ----------
    distances : np.ndarray
        array containing optical distances values
    data : np.ndarray
        array containing data to plot.
    y_label : str, optional
        y-axis label. The default is ''.
    x_label : str, optional
        x-axis label. The default is ''.
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
    win.setWindowIcon(QtGui.QIcon(str(Path(__file__).resolve().parent / 'logo.png')))
    win.setBackground('w') # white bg
    win.resize(1200, 500)

    plot = win.addPlot(title=title)
    plot.setTitle(title, size='20pt', color='k')
    plot.setLabel('left', y_label, **{'color': 'k', 'font-size': '14pt'})
    plot.plot(distances, data, pen=pg.mkPen('k', width=1))

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
    plot.setXRange(min(distances), max(distances), padding=0)
    plot.getViewBox().setLimits(xMin=min(distances), xMax=max(distances),
                                yMin=min(data), yMax=max(data))

    pg.exec()



'''
-----------------------------------------------------------------
Matrix Plot Functions
-----------------------------------------------------------------
'''

def plot_2d_timeseries(timestamps: np.ndarray, data: np.ndarray, y_ticks: list,
                       max_value: float = None, y_label: str = '',
                       title: str = '', cmap:str = 'seismic',
                       cbar_label:str = '') -> None :
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

    Returns
    -------
    None
        -.
    '''

    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication(sys.argv)

    win = pg.GraphicsLayoutWidget(show=True)
    win.setWindowTitle('Fobench Data Plot')
    win.setWindowIcon(QtGui.QIcon(str(Path(__file__).resolve().parent / 'logo.png')))
    win.setBackground('w')
    win.resize(1200, 800)

    plot = win.addPlot(title=title)
    plot.setTitle(title, size='20pt', color='k')
    plot.setLabel('left', y_label, **{'color': 'k', 'font-size': '14pt'})

    plot.setCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))
    plot.setMouseEnabled(x=True, y=True)
    plot.getViewBox().setMouseMode(pg.ViewBox.RectMode) # One-button

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
    plot.setAxisItems({'bottom': pg.DateAxisItem(utcOffset=1)}) # when is utcoffset necessary??
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

    # cursor tracking in data units
    label = pg.LabelItem(justify='left')
    win.addItem(label, row=2, col=0)

    def mouse_moved(evt):
        pos = evt[0]
        if plot.sceneBoundingRect().contains(pos):
            mouse_point = plot.vb.mapSceneToView(pos)
            x_val = mouse_point.x()
            y_val = mouse_point.y()
            x_datetime = datetime.datetime.utcfromtimestamp(x_val).strftime('%Y-%m-%d %H:%M:%S')
            label.setText(f'Time: {x_datetime} | {y_label}: {y_val:.1f}', color='k')
    proxy = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60, slot=mouse_moved)

    values = (-max_value, max_value) if max_value is not None else None
    data_range = np.nanmax(data) - np.nanmin(data)

    cmap = pg.colormap.get(cmap, source='matplotlib')
    bar = pg.ColorBarItem(colorMap=cmap, values=values,
                          label=cbar_label, interactive=True, rounding=0.001*data_range)
    bar.setImageItem(image, insert_in=plot)

    pg.exec()


def plot_2d_distance(distances: np.ndarray, data: np.ndarray, y_ticks: list,
                     max_value: float = None, y_label: str = '',
                     x_label: str = 'Optical Distance [m]',
                     title: str = '', cmap:str = 'seismic',
                     cbar_label:str = '', invert_y=False) -> None :
    '''
    generate generic matrix plot where x-axis represents distance

    Parameters
    ----------
    distances : np.ndarray
        array containing optical distances of channels.
    data : np.ndarray
        array containing data to plot.
    y_ticks : list
        y-axis tick labels.
    max_value : float, optional
        value that sets limits of colorbar. The default is None.
    y_label : str, optional
        y-axis label. The default is ''.
    x_label : str, optional
        x-axis label. The default is 'Optical Distance [m]'.
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
    win.setWindowIcon(QtGui.QIcon(str(Path(__file__).resolve().parent / 'logo.png')))
    win.setBackground('w')
    win.resize(1200, 800)

    plot = win.addPlot(title=title)
    plot.setTitle(title, size='20pt', color='k')
    plot.setLabel('left', y_label, **{'color': 'k', 'font-size': '14pt'})
    plot.setLabel('bottom', x_label, **{'color': 'k', 'font-size': '14pt'})

    plot.setCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))
    plot.setMouseEnabled(x=True, y=True)
    plot.getViewBox().setMouseMode(pg.ViewBox.RectMode) # One-button
    if invert_y:
        plot.getViewBox().invertY(True)

    x_min, x_max = distances[0], distances[-1]
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

    plot.setXRange(min(distances), max(distances), padding=0)
    plot.setYRange(min(y_ticks), max(y_ticks), padding=0)
    plot.getViewBox().setLimits(xMin=min(distances), xMax=max(distances),
                                yMin=min(y_ticks), yMax=max(y_ticks))

    # cursor tracking in data units
    label = pg.LabelItem(justify='left')
    win.addItem(label, row=2, col=0)

    def mouse_moved(evt):
        pos = evt[0]
        if plot.sceneBoundingRect().contains(pos):
            mouse_point = plot.vb.mapSceneToView(pos)
            x_val = mouse_point.x()
            y_val = mouse_point.y()
            label.setText(f'{x_label}: {x_val:.1f} | {y_label}: {y_val:.1f}', color='k')
    proxy = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60, slot=mouse_moved)

    values = (-max_value, max_value) if max_value is not None else None
    data_range = np.nanmax(data) - np.nanmin(data)

    cmap = pg.colormap.get(cmap, source='matplotlib')
    bar = pg.ColorBarItem(colorMap=cmap, values=values,
                          label=cbar_label, interactive=True, rounding=0.001*data_range)
    bar.setImageItem(image, insert_in=plot)

    pg.exec()