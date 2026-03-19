'''
Contains all functionality related to plotting using PyQtGraph, i.e.
whenever plot_mode is set to 'pyqt'
'''
import sys
import datetime
from pathlib import Path
import numpy as np
from PyQt5 import QtCore
import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets, QtGui

'''Line Plot Functions'''

def plot_timeseries(timestamps: np.ndarray, data: np.ndarray, dt: float,
                    y_label: str = '', title: str = '') -> None:
    '''
    Generate generic time series plot using PyQtGraph, ideal for channel plots

    Parameters
    ----------
    timestamps : np.ndarray
        Array containing Unix timestamps of data.
    data : np.ndarray
        Array containing data to plot.
    dt : float
        Sampling period of data.
    y_label : str, optional
        y-axis label.
    title : str, optional
        Title of plot..

    Returns
    -------
    None
        -
    '''

    win, app, plot, y_axis, x_axis = get_layout(size=(1200, 500), win_title=title,
                                                x_is_time=True)
    plot.setLabel('left', y_label, **{'color': 'k', 'font-size': '14pt'})
    plot.plot(timestamps, data, pen=pg.mkPen('k', width=1))
    date = datetime.datetime.fromtimestamp(timestamps[0]).strftime('%d.%m.%Y')
    x_axis.setLabel(date, **{'color': 'k', 'font-size': '14pt'})
    plot.setXRange(min(timestamps), max(timestamps), padding=0)
    plot.getViewBox().setLimits(xMin=min(timestamps), xMax=max(timestamps),
                                yMin=min(data), yMax=max(data))
    label = pg.LabelItem(justify='left', size='10pt', color='black')
    win.addItem(label, row=2, col=0)
    mouse_moved = tracker_factory(plot=plot, label=label, dt=dt, label_text='Time: {x}')
    proxy = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60, slot=mouse_moved)
    pg.exec()

def plot_record_section(timestamps: np.ndarray, data: np.ndarray, dt: float,
                        numbers: np.ndarray, y_label: str = '', title: str = '') -> None:
    '''
    Extended version of timeseries plot for multi-channel data

    Parameters
    ----------
    timestamps : np.ndarray
        Array containing Unix timestamps of data.
    data : np.ndarray
        Array containing data to plot.
    dt : float
        Sampling period of data.
    numbers : np.ndarray
        Array containing channel numbers for plotting.
    y_label : str, optional
        y-axis label.
    title : str, optional
        Title of plot.

    Returns
    -------
    None
        -
    '''
    win, app, plot, y_axis, x_axis = get_layout(size=(1200, 500), win_title=title,
                                                x_is_time=True)
    plot.setLabel('left', y_label, **{'color': 'k', 'font-size': '14pt'})
    data, numbers = data[:, ::-1], numbers[::-1]
    offset = data.max()
    ticks = []
    for i in range(data.shape[1]):
        y = -i*offset
        plot.plot(timestamps, data[:, i]+y, pen=pg.mkPen('k', width=1))
        ticks.append((y, str(numbers[i])))
    y_axis.setTicks([ticks])
    y_axis.enableAutoSIPrefix(False)
    date = datetime.datetime.fromtimestamp(timestamps[0]).strftime('%d.%m.%Y')
    x_axis.setLabel(date, **{'color': 'k', 'font-size': '14pt'})
    plot.setXRange(min(timestamps), max(timestamps), padding=0)
    plot.getViewBox().setLimits(xMin=min(timestamps), xMax=max(timestamps),
                                yMax=data[:, 0].max(),
                                yMin=data[:, -1].min()-(data.shape[1]-1)*offset)
    label = pg.LabelItem(justify='left', size='10pt', color='black')
    win.addItem(label, row=2, col=0)
    mouse_moved = tracker_factory(plot=plot, label=label, dt=dt, label_text='Time: {x}')
    proxy = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60, slot=mouse_moved)
    pg.exec()


def plot_distance(distances: np.ndarray, channels_num: np.ndarray, data: np.ndarray,
                  y_label: str = '', x_label: str = 'Channel',
                  title: str = '') -> None:
    '''
    Generate generic distances series plot

    Parameters
    ----------
    distances : np.ndarray
        Array containing optical distances values.
    channels_num : np.ndarray
        Array containing channel numbers.
    data : np.ndarray
        Array containing data to plot.
    y_label : str, optional
        y-axis label.
    x_label : str, optional
        x-axis label.
    title : str, optional
        Title of plot.

    Returns
    -------
    None
        -
    '''
    win, app, plot, y_axis, x_axis = get_layout(size=(1200, 500), win_title=title)
    plot.setLabel('left', y_label, **{'color': 'k', 'font-size': '14pt'})
    plot.setLabel('bottom', x_label, **{'color': 'k', 'font-size': '14pt'})
    text_label, h_layout, container = get_bottom_layout()
    h_layout.addWidget(button := get_axis_button())
    state = {'x_vals': channels_num, 'dx': channels_num[1]-channels_num[0]}

    def refresh(x_vals):
        state['x_vals'], state['dx'] = x_vals, x_vals[1]-x_vals[0]
        dx = state['dx']
        plot.clear()
        plot.plot(x_vals, data, pen=pg.mkPen('k', width=1))
        plot.getViewBox().setLimits(xMin=-np.inf, xMax=np.inf)
        plot.setXRange(min(x_vals), max(x_vals), padding=0)
        plot.getViewBox().setLimits(xMin=min(x_vals), xMax=max(x_vals))
        state['proxy'] = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60,
                                        slot=tracker_factory(plot=plot, label=text_label, dx=dx,
                                                             label_text=x_axis.labelText+': {x}'))
    def switch_axis():
        if button.text() == 'Distance':
            button.setText('Channel')
            plot.setLabel('bottom', 'Optical Distance [m]', **{'color': 'k', 'font-size': '14pt'})
            refresh(distances)
        else:
            button.setText('Distance')
            plot.setLabel('bottom', 'Channel', **{'color': 'k', 'font-size': '14pt'})
            refresh(channels_num)

    button.clicked.connect(switch_axis)
    refresh(channels_num)
    proxy_container = QtWidgets.QGraphicsProxyWidget()
    proxy_container.setWidget(container)
    win.addItem(proxy_container, row=2, col=0)
    pg.exec()

def plot_spectral(frequencies: np.ndarray, amplitudes: np.ndarray,
                  y_label: str = 'Amplitude', x_label: str = 'Frequency [Hz]',
                  title: str = '') -> None:
    '''
    Generate generic amplitude over frequency plot

    Parameters
    ----------
    frequencies : np.ndarray
        Array containing frequency values.
    amplitudes : np.ndarray
        Array containing amplitudes to plot.
    y_label : str, optional
        y-axis label.
    x_label : str, optional
        x-axis label.
    title : str, optional
        Title of plot.

    Returns
    -------
    None
        -
    '''

    win, app, plot, y_axis, x_axis = get_layout(size=(1200, 500), win_title=title)
    plot.setLabel('left', y_label, **{'color': 'k', 'font-size': '14pt'})
    plot.setLabel('bottom', x_label, **{'color': 'k', 'font-size': '14pt'})
    dx = frequencies[1]-frequencies[0]
    plot.plot(frequencies, amplitudes, pen=pg.mkPen('k', width=1))
    plot.setXRange(min(frequencies), max(frequencies), padding=0)
    plot.getViewBox().setLimits(xMin=min(frequencies), xMax=max(frequencies),
                                yMin=min(amplitudes), yMax=max(amplitudes))
    text_label, h_layout, container = get_bottom_layout()
    mouse_moved = tracker_factory(plot=plot, label=text_label, dx=dx, label_text=x_axis.labelText+': {x:.2f}')
    proxy = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60, slot=mouse_moved)
    proxy_container = QtWidgets.QGraphicsProxyWidget()
    proxy_container.setWidget(container)
    win.addItem(proxy_container, row=2, col=0)
    pg.exec()

'''Matrix Plot Functions'''

def plot_2d_timeseries(timestamps: np.ndarray, data: np.ndarray, y_ticks: list,
                       dt: float, vmin: float = None, vmax: float = None, y_label: str = '',
                       title: str = '', cmap: str = 'seismic',
                       cbar_label: str = '', distances: np.ndarray = None) -> None:
    '''
    Generate generic matrix plot where x-axis represents time, e.g.
    waterfall visualisation of data or spectrograms, PSDs, RMSA ...

    Parameters
    ----------
    timestamps : np.ndarray
        Array containing Unix timestamps of data.
    data : np.ndarray
        Array containing data to plot.
    y_ticks : list
        y-axis tick labels.
    dt : float
        Sampling period of data.
    vmin, vmax : float, optional
        Minimum and maximum limits of colorbar.
    y_label : str, optional
        y-axis label.
    title : str, optional
        Plot title.
    cmap : str, optional
        Colormap to use. The default is 'seismic'.
    cbar_label : str, optional
        Label of colorbar.
    distances : np.ndarray, optional
        Optional distances to use for y-axis.

    Returns
    -------
    None
        -.
    '''
    win, app, plot, y_axis, x_axis = get_layout(size=(1200, 800), win_title=title,
                                                 x_is_time=True)
    plot.setLabel('left', y_label, **{'color': 'k', 'font-size': '14pt'})
    plot.setCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))
    plot.setMouseEnabled(x=True, y=True)
    plot.getViewBox().setMouseMode(pg.ViewBox.RectMode)
    x_min, x_max = timestamps[0], timestamps[-1]
    image = pg.ImageItem()
    plot.addItem(image)
    image.setImage(data)
    text_label, h_layout, container = get_bottom_layout()
    h_layout.addWidget(scale_button := get_colorscale_button())
    state = {'y_vals':y_ticks, 'dy':y_ticks[1]-y_ticks[0]}

    def refresh(y_vals):
        state['y_vals'], state['dy'] = y_vals, y_vals[1]-y_vals[0]
        dy = state['dy']
        y_min, y_max = y_vals[0], y_vals[-1]
        image.setRect(x_min, y_min-dy/2, x_max-x_min, y_max-y_min+dy)
        plot.setXRange(x_min, x_max, padding=0)
        plot.getViewBox().setLimits(xMin=x_min, xMax=x_max,
                                    yMin=-np.inf, yMax=np.inf)
        plot.setYRange(min(y_vals)-dy/2, max(y_vals)+dy/2, padding=0)
        plot.getViewBox().setLimits(xMin=x_min, xMax=x_max,
                                    yMin=min(y_vals)-dy/2, yMax=max(y_vals)+dy/2)
        state['proxy'] = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60,
                                        slot=tracker_factory(plot=plot, label=text_label, dt=dt, dy=dy,
                                                             label_text='Time: {x} | '+y_axis.labelText+': {y}'))
    if distances is not None:
        h_layout.addWidget(button := get_axis_button())

        def switch_axis():
            if button.text() == 'Distance':
                button.setText('Channel')
                plot.setLabel('left', 'Optical Distance [m]', **{'color': 'k', 'font-size': '14pt'})
                refresh(distances)
            else:
                button.setText('Distance')
                plot.setLabel('left', 'Channel', **{'color': 'k', 'font-size': '14pt'})
                refresh(y_ticks)

        button.clicked.connect(switch_axis)
    refresh(y_ticks)
    data_range = np.nanmax(data)-np.nanmin(data)
    cmap = pg.colormap.get(cmap, source='matplotlib')
    bar = pg.ColorBarItem(colorMap=cmap, values=(vmin, vmax),
                          label=cbar_label, interactive=True, rounding=0.0001*data_range)
    bar.setImageItem(image, insert_in=plot)
    scale_button.clicked.connect(lambda: reset_scale(bar, vmin, vmax))
    proxy_container = QtWidgets.QGraphicsProxyWidget()
    proxy_container.setWidget(container)
    win.addItem(proxy_container, row=2, col=0)
    pg.exec()


def plot_2d_distance(distances: np.ndarray, channels_num: np.ndarray,
                     data: np.ndarray, y_ticks: list, vmin: float = None, vmax: float = None,
                     y_label: str = '', x_label: str = 'Channel', title: str = '',
                     cmap: str = 'seismic', cbar_label: str = '', invert_y=False) -> None:
    '''
    Generate generic matrix plot where x-axis represents distance.

    Parameters
    ----------
    distances : np.ndarray
        Array containing optical distances of channels.
    channels_num : np.ndarray
        Array containing channel numbers.
    data : np.ndarray
        Array containing data to plot.
    y_ticks : list
        y-axis tick labels.
    vmin, vmax : float, optional
        Minimum and maximum limits of colorbar. The default is None.
    y_label : str, optional
        y-axis label. The default is ''.
    x_label : str, optional
        x-axis label. The default is 'Channel'.
    title : str, optional
        Title of plot. The default is ''.
    cmap : str, optional
        Colormap to use. The default is 'seismic'.
    cbar_label : str, optional
        Label of colorbar. The default is ''.
    invert_y : bool, optional
        Invert y-axis. The default is False.

    Returns
    -------
    None
        -.
    '''
    win, app, plot, y_axis, x_axis = get_layout(size=(1200, 800), win_title=title)
    dy = y_ticks[1]-y_ticks[0]
    y_min, y_max = y_ticks[0], y_ticks[-1]
    plot.setLabel('left', y_label, **{'color': 'k', 'font-size': '14pt'})
    plot.setLabel('bottom', x_label, **{'color': 'k', 'font-size': '14pt'})
    plot.setCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))
    plot.setMouseEnabled(x=True, y=True)
    plot.getViewBox().setMouseMode(pg.ViewBox.RectMode)
    if invert_y:
        plot.getViewBox().invertY(True)
    image = pg.ImageItem()
    plot.addItem(image)
    image.setImage(data)
    plot.setYRange(min(y_ticks)-dy/2, max(y_ticks)+dy/2, padding=0)
    plot.getViewBox().setLimits(yMin=min(y_ticks)-dy/2, yMax=max(y_ticks)+dy/2)
    text_label, h_layout, container = get_bottom_layout()
    h_layout.addWidget(button := get_axis_button())
    h_layout.addWidget(scale_button := get_colorscale_button())
    state = {'x_vals': channels_num, 'dx': channels_num[1]-channels_num[0]}

    def refresh(x_vals):
        state['x_vals'], state['dx'] = x_vals, x_vals[1]-x_vals[0]
        dx = state['dx']
        x_min, x_max = x_vals[0], x_vals[-1]
        image.setRect(x_min-dx/2, y_min-dy/2, x_max-x_min+dx, y_max-y_min+dy)
        plot.getViewBox().setLimits(xMin=-np.inf, xMax=np.inf,
                                    yMin=min(y_ticks)-dy/2, yMax=max(y_ticks)+dy/2)
        plot.setXRange(min(x_vals)-dx/2, max(x_vals)+dx/2, padding=0)
        plot.getViewBox().setLimits(xMin=min(x_vals)-dx/2, xMax=max(x_vals)+dx/2,
                                    yMin=min(y_ticks)-dy/2, yMax=max(y_ticks)+dy/2)
        state['proxy'] = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60,
                                        slot=tracker_factory(plot=plot, label=text_label, dx=dx, dy=dy,
                                                             label_text=x_axis.labelText+': {x} | '+y_label+': {y:.2f}'))

    def switch_axis():
        if button.text() == 'Distance':
            button.setText('Channel')
            plot.setLabel('bottom', 'Optical Distance [m]', **{'color': 'k', 'font-size': '14pt'})
            refresh(distances)
        else:
            button.setText('Distance')
            plot.setLabel('bottom', 'Channel', **{'color': 'k', 'font-size': '14pt'})
            refresh(channels_num)

    button.clicked.connect(switch_axis)
    refresh(channels_num)
    data_range = np.nanmax(data)-np.nanmin(data)
    cmap = pg.colormap.get(cmap, source='matplotlib')
    bar = pg.ColorBarItem(colorMap=cmap, values=(vmin, vmax), label=cbar_label,
                          interactive=True, rounding=0.0001*data_range)
    bar.setImageItem(image, insert_in=plot)

    scale_button.clicked.connect(lambda: reset_scale(bar, vmin, vmax))

    proxy_container = QtWidgets.QGraphicsProxyWidget()
    proxy_container.setWidget(container)
    win.addItem(proxy_container, row=2, col=0)
    pg.exec()

'''Helper Functions'''

def get_layout(size: tuple = (1200, 600), win_title: str = None,
               x_is_time: bool = False):
    '''Returns common basic layout of plots'''
    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication(sys.argv)

    win = pg.GraphicsLayoutWidget(show=True)
    win.setWindowTitle(f'Fobench: {win_title}')
    win.setWindowIcon(QtGui.QIcon(
        str(Path(__file__).resolve().parent / 'logo.png')))
    win.setBackground('w')
    win.resize(*size)
    plot = win.addPlot(title=win_title)
    plot.setTitle(win_title, size='18pt', color='k')

    if x_is_time:
        plot.setAxisItems({'bottom': pg.DateAxisItem(utcOffset=1)})
    for axis in (y_axis := plot.getAxis('left'), x_axis := plot.getAxis('bottom')):
        axis.setPen(pg.mkPen('k', width=2))
        axis.setTextPen(pg.mkPen('k'))
        axis.setStyle(tickFont=pg.Qt.QtGui.QFont('Arial', 14))
    return win, app, plot, y_axis, x_axis

def tracker_factory(plot: pg.PlotItem, label: pg.LabelItem, label_text: str,
                    dt:float = None, dx:float = None, dy: float = None):
    '''Factory that returns a function tracking the mouse position and displays it in data units'''
    x_step = dt or dx
    def mouse_moved(evt):
        pos = evt[0]
        if plot.sceneBoundingRect().contains(pos):
            mp = plot.vb.mapSceneToView(pos)
            x = round(mp.x()/x_step)*x_step
            y = round(mp.y()/dy)*dy if dy is not None else mp.y()
            if dt is not None:
                x = datetime.datetime.utcfromtimestamp(x).strftime(
                    '%Y-%m-%d %H:%M:%S.%f')
            label.setText(label_text.format(x=x, y=y))
    return mouse_moved

def get_axis_button() -> QtWidgets.QPushButton:
    '''Returns a button that can be used for switching axis between optical distance and channel number'''
    button = QtWidgets.QPushButton('Distance')
    button.setToolTip('Switch x Axis between Channel Number and Optical Distance')
    button.setFixedWidth(70)
    return button

def get_colorscale_button() -> QtWidgets.QPushButton:
    '''Returns a button that can be used for resetting colorscale'''
    button = QtWidgets.QPushButton('Reset colorscale')
    button.setToolTip('Resets the colorscale to initial values')
    button.setFixedWidth(110)
    return button

def reset_scale(bar, vmin, vmax) -> None:
    '''Function to resest the colorscale to original values'''
    bar.setLevels(values=(vmin, vmax))

def get_bottom_layout():
    '''Returns bottom layout that holds label for tracking mouse'''
    container = QtWidgets.QWidget()
    h_layout = QtWidgets.QHBoxLayout(container)
    h_layout.setContentsMargins(5, 5, 5, 5)
    text_label = QtWidgets.QLabel()
    text_label.setStyleSheet('color: black; font-size: 10pt;')
    h_layout.addWidget(text_label)
    return text_label, h_layout, container