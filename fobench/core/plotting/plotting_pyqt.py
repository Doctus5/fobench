"""Contains all functionality related to plotting using PyQtGraph, i.e.
whenever plot_mode is set to 'pyqt'."""

import sys
import datetime
from pathlib import Path
import numpy as np
from PyQt5 import QtCore
import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets, QtGui
import matplotlib.pyplot as plt

"""Line Plot Functions"""

def plot_timeseries(timestamps: np.ndarray, data: np.ndarray | list, dt: float,
                    y_label: str = "", title: str = "", labels: list[str] | None = None) -> None:

    """Generate generic time series plot using PyQtGraph, ideal for channel plots.

    Parameters
    ----------
    timestamps : np.ndarray
        Array containing Unix timestamps of data.
    data : np.ndarray | list
        Array containing data to plot.
    dt : float
        Sampling period of data.
    y_label : str, optional
        y-axis label.
    title : str, optional
        Title of plot.
    labels : list[str] | None
        Legend labels for each amplitude array. If None, no legend is shown.

    Returns
    -------
    None

    """

    y_label = y_label.title()
    win, app, plot, y_axis, x_axis = get_layout(size=(1200, 500), win_title=title,
                                                x_is_time=True)
    plot.setLabel("left", y_label, **{"color": "k", "font-size": "14pt"})
    if isinstance(data, np.ndarray) and data.ndim == 1:
        data = [data]
    elif isinstance(data, np.ndarray) and data.ndim == 2:
        data = [data[:, i] for i in range(data.shape[1])]
    if labels is not None:
        plot.addLegend(pen="k", labelTextColor="k")
    colors = ["k"] if len(data) == 1 else get_colors(len(data), colormap="tab10")
    for i, channel in enumerate(data):
        label = labels[i] if labels is not None else None
        plot.plot(timestamps, channel, pen=pg.mkPen(colors[i], width=1), name=label)
    all_data = np.concatenate(data)
    date = datetime.datetime.fromtimestamp(timestamps[0]).strftime("%d.%m.%Y")
    x_axis.setLabel(date, **{"color": "k", "font-size": "14pt"})
    plot.setXRange(min(timestamps), max(timestamps), padding=0)
    plot.getViewBox().setLimits(xMin=min(timestamps), xMax=max(timestamps),
                                yMin=min(all_data), yMax=max(all_data))
    label_item = pg.LabelItem(justify="left", size="10pt", color="black")
    win.addItem(label_item, row=2, col=0)
    label_text = "Time: {x} | "+y_label+": {y:1e}"
    mouse_moved = tracker_factory(plot=plot, label=label_item, dt=dt, label_text=label_text)
    proxy = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60, slot=mouse_moved)
    pg.exec()

def plot_record_section(timestamps: np.ndarray, data: np.ndarray, dt: float,
                        numbers: np.ndarray, y_label: str = "", title: str = "") -> None:

    """Extended version of timeseries plot for multi-channel data.

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

    """

    win, app, plot, y_axis, x_axis = get_layout(size=(1200, 500), win_title=title,
                                                x_is_time=True)
    plot.setLabel("left", y_label, **{"color": "k", "font-size": "14pt"})
    data, numbers = data[:, ::-1], numbers[::-1]
    offset = data.max()
    ticks = []
    for i in range(data.shape[1]):
        y = -i*offset
        plot.plot(timestamps, data[:, i]+y, pen=pg.mkPen("k", width=1))
        ticks.append((y, str(numbers[i])))
    y_axis.setTicks([ticks])
    y_axis.enableAutoSIPrefix(False)
    date = datetime.datetime.fromtimestamp(timestamps[0]).strftime("%d.%m.%Y")
    x_axis.setLabel(date, **{"color": "k", "font-size": "14pt"})
    plot.setXRange(min(timestamps), max(timestamps), padding=0)
    plot.getViewBox().setLimits(xMin=min(timestamps), xMax=max(timestamps),
                                yMax=data[:, 0].max(),
                                yMin=data[:, -1].min()-(data.shape[1]-1)*offset)
    label = pg.LabelItem(justify="left", size="10pt", color="black")
    win.addItem(label, row=2, col=0)
    mouse_moved = tracker_factory(plot=plot, label=label, dt=dt, label_text="Time: {x}")
    proxy = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60, slot=mouse_moved)
    pg.exec()


def plot_distance(distances: np.ndarray, channels_num: np.ndarray, data: np.ndarray,
                  y_label: str = "", x_label: str = "Channel",
                  title: str = "") -> None:

    """Generate generic distances series plot.

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

    """

    win, app, plot, y_axis, x_axis = get_layout(size=(1200, 500), win_title=title)
    plot.setLabel("left", y_label, **{"color": "k", "font-size": "14pt"})
    plot.setLabel("bottom", x_label, **{"color": "k", "font-size": "14pt"})
    text_label, h_layout, container = get_bottom_layout()
    h_layout.addWidget(button := get_axis_button())
    state = {"x_vals": channels_num, "dx": channels_num[1]-channels_num[0]}
    def refresh(x_vals):
        state["x_vals"], state["dx"] = x_vals, x_vals[1]-x_vals[0]
        dx = state["dx"]
        plot.clear()
        plot.plot(x_vals, data, pen=pg.mkPen("k", width=1))
        plot.getViewBox().setLimits(xMin=-np.inf, xMax=np.inf)
        plot.setXRange(min(x_vals), max(x_vals), padding=0)
        plot.getViewBox().setLimits(xMin=min(x_vals), xMax=max(x_vals),
                                    yMin=min(data), yMax=max(data))
        label_text = x_axis.labelText+": {x} | "+y_label+": {y:1e}"
        state["proxy"] = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60,
                                        slot=tracker_factory(plot=plot, label=text_label, dx=dx,
                                                             label_text=label_text))
    def switch_axis():
        if button.text() == "Distance":
            button.setText("Channel")
            plot.setLabel("bottom", "Optical Distance [m]", **{"color": "k", "font-size": "14pt"})
            refresh(distances)
        else:
            button.setText("Distance")
            plot.setLabel("bottom", "Channel", **{"color": "k", "font-size": "14pt"})
            refresh(channels_num)

    button.clicked.connect(switch_axis)
    refresh(channels_num)
    proxy_container = QtWidgets.QGraphicsProxyWidget()
    proxy_container.setWidget(container)
    win.addItem(proxy_container, row=2, col=0)
    pg.exec()

def plot_spectral(frequencies: np.ndarray, amplitudes: list[np.ndarray] | np.ndarray,
                  y_label: str = "Amplitude", x_label: str = "Frequency [Hz]",
                  title: str = "", labels: list[str] | None = None) -> None:

    """Generate generic amplitude over frequency plot.

    Parameters
    ----------
    frequencies : np.ndarray
        Array containing frequency values.
    amplitudes : list[np.ndarray] | np.ndarray
        Single array or list of arrays containing amplitudes to plot.
    y_label : str, optional
        y-axis label.
    x_label : str, optional
        x-axis label.
    title : str, optional
        Title of plot.
    labels : list[str] | None, optional
        Legend labels for each amplitude array. If None, no legend is shown.
    Returns
    -------
    None

    """

    win, app, plot, y_axis, x_axis = get_layout(size=(1200, 500), win_title=title)
    plot.setLabel("left", y_label, **{"color": "k", "font-size": "14pt"})
    plot.setLabel("bottom", x_label, **{"color": "k", "font-size": "14pt"})
    dx = frequencies[1] - frequencies[0]
    if labels is not None:
        plot.addLegend(pen="k", labelTextColor="k")
    all_amplitudes = np.concatenate(amplitudes)
    n = amplitudes.shape[1] if isinstance(amplitudes, np.ndarray) else len(amplitudes)
    colors = ["k"] if n == 1 else get_colors(n, colormap="tab10")
    for i, amp in enumerate(amplitudes.T):
        label = labels[i] if labels is not None else None
        plot.plot(frequencies, amp, pen=pg.mkPen(colors[i], width=1), name=label)
    plot.setXRange(min(frequencies), max(frequencies), padding=0)
    plot.getViewBox().setLimits(xMin=min(frequencies), xMax=max(frequencies),
                                yMin=min(all_amplitudes), yMax=max(all_amplitudes))
    text_label, h_layout, container = get_bottom_layout()
    label_text = x_axis.labelText + ": {x:.2f} | "+y_label+": {y:e}"
    mouse_moved = tracker_factory(plot=plot, label=text_label, dx=dx, label_text=label_text)
    proxy = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60, slot=mouse_moved)
    proxy_container = QtWidgets.QGraphicsProxyWidget()
    proxy_container.setWidget(container)
    win.addItem(proxy_container, row=2, col=0)
    pg.exec()


"""Matrix Plot Functions"""

def plot_2d_timeseries(timestamps: np.ndarray, data: np.ndarray, y_ticks: list,
                       dt: float, vmin: float = None, vmax: float = None, y_label: str = "",
                       title: str = "", cmap: str = "seismic",
                       cbar_label: str = "", distances: np.ndarray = None) -> None:

    """Generate generic matrix plot where x-axis represents time.

    Parameters
    ----------
    timestamps : np.ndarray
        Array containing timestamps.
    data : np.ndarray
        Data to plot.
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
        Colormap to be used.
    cbar_label : str, optional
        Colorbar label
    distances : np.ndarray, optional
        Array containing optical distance values.

    Returns
    -------
    None

    """

    cbar_label = cbar_label.title()

    win, app, plot, y_axis, x_axis = get_layout(size=(1200, 800), win_title=title,
                                                 x_is_time=True)
    plot.setLabel("left", y_label, **{"color": "k", "font-size": "14pt"})
    plot.setCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))
    plot.setMouseEnabled(x=True, y=True)
    plot.getViewBox().setMouseMode(pg.ViewBox.RectMode)
    x_min, x_max = timestamps[0], timestamps[-1]
    image = pg.ImageItem()
    plot.addItem(image)
    image.setImage(data)
    text_label, h_layout, container = get_bottom_layout()
    state = {"y_vals": y_ticks, "dy": y_ticks[1] - y_ticks[0]}

    def refresh(y_vals):
        state["y_vals"], state["dy"] = y_vals, y_vals[1] - y_vals[0]
        dy = state["dy"]
        y_min, y_max = y_vals[0], y_vals[-1]
        image.setRect(x_min, y_min - dy / 2, x_max - x_min, y_max - y_min + dy)
        plot.setXRange(x_min, x_max, padding=0)
        plot.getViewBox().setLimits(xMin=x_min, xMax=x_max, yMin=-np.inf, yMax=np.inf)
        plot.setYRange(min(y_vals) - dy / 2, max(y_vals) + dy / 2, padding=0)
        plot.getViewBox().setLimits(xMin=x_min, xMax=x_max,
                                    yMin=min(y_vals) - dy / 2, yMax=max(y_vals) + dy / 2)
        state["proxy"] = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60,
                                        slot=tracker_factory(plot=plot, label=text_label, dt=dt, dy=dy,
                                                             label_text="Time: {x} | " + y_axis.labelText + ": {y} | " + cbar_label + " : {z:2e}",
                                                             item=image))

    if distances is not None:
        h_layout.addWidget(button := get_axis_button())
        def switch_axis():
            if button.text() == "Distance":
                button.setText("Channel")
                plot.setLabel("left", "Optical Distance [m]", **{"color": "k", "font-size": "14pt"})
                refresh(distances)
            else:
                button.setText("Distance")
                plot.setLabel("left", "Channel", **{"color": "k", "font-size": "14pt"})
                refresh(y_ticks)
        button.clicked.connect(switch_axis)
    refresh(y_ticks)
    data_range = np.nanmax(data) - np.nanmin(data)
    data_min, data_max = np.nanmin(data), np.nanmax(data)

    cmap_obj = pg.colormap.get(cmap, source='matplotlib')
    _vmin = vmin if vmin is not None else data_min
    _vmax = vmax if vmax is not None else data_max
    bar = pg.ColorBarItem(colorMap=cmap_obj, values=(_vmin, _vmax),
                          label=cbar_label, interactive=True, rounding=0.0001 * data_range)
    bar.setImageItem(image, insert_in=plot)
    h_layout.addWidget(scale_button := get_colorscale_button())
    scale_button.clicked.connect(lambda: reset_scale(bar, vmin, vmax))
    h_layout.addWidget(get_vminmax_button(bar))

    proxy_container = QtWidgets.QGraphicsProxyWidget()
    proxy_container.setWidget(container)
    win.addItem(proxy_container, row=2, col=0)
    pg.exec()

def plot_2d_distance(distances: np.ndarray, channels_num: np.ndarray,
                     data: np.ndarray, y_ticks: list, vmin: float = None, vmax: float = None,
                     y_label: str = "", x_label: str = "Channel", title: str = "",
                     cmap: str = "seismic", cbar_label: str = "", invert_y=False) -> None:

    """Generate generic matrix plot where x-axis represents distance.

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
        Minimum and maximum limits of colorbar.
    y_label : str, optional
        y-axis label.
    x_label : str, optional
        x-axis label.
    title : str, optional
        Title of plot.
    cmap : str, optional
        Colormap to use.
    cbar_label : str, optional
        Label of colorbar.
    invert_y : bool, optional
        Invert y-axis.

    Returns
    -------
    None

    """

    cbar_label=cbar_label.title()
    win, app, plot, y_axis, x_axis = get_layout(size=(1200, 800), win_title=title)
    dy = y_ticks[1]-y_ticks[0]
    y_min, y_max = y_ticks[0], y_ticks[-1]
    plot.setLabel("left", y_label, **{"color": "k", "font-size": "14pt"})
    plot.setLabel("bottom", x_label, **{"color": "k", "font-size": "14pt"})
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
    state = {"x_vals": channels_num, "dx": channels_num[1]-channels_num[0]}

    def refresh(x_vals):
        state["x_vals"], state["dx"] = x_vals, x_vals[1]-x_vals[0]
        dx = state["dx"]
        x_min, x_max = x_vals[0], x_vals[-1]
        image.setRect(x_min-dx/2, y_min-dy/2, x_max-x_min+dx, y_max-y_min+dy)
        plot.getViewBox().setLimits(xMin=-np.inf, xMax=np.inf,
                                    yMin=min(y_ticks)-dy/2, yMax=max(y_ticks)+dy/2)
        plot.setXRange(min(x_vals)-dx/2, max(x_vals)+dx/2, padding=0)
        plot.getViewBox().setLimits(xMin=min(x_vals)-dx/2, xMax=max(x_vals)+dx/2,
                                    yMin=min(y_ticks)-dy/2, yMax=max(y_ticks)+dy/2)
        label_text = x_axis.labelText+": {x} | "+y_label+": {y:.2f} | "+cbar_label+": {z:e}"
        state["proxy"] = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60,
                                        slot=tracker_factory(plot=plot, label=text_label, dx=dx, dy=dy,
                                                             label_text=label_text, item=image))

    def switch_axis():
        if button.text() == "Distance":
            button.setText("Channel")
            plot.setLabel("bottom", "Optical Distance [m]", **{"color": "k", "font-size": "14pt"})
            refresh(distances)
        else:
            button.setText("Distance")
            plot.setLabel("bottom", "Channel", **{"color": "k", "font-size": "14pt"})
            refresh(channels_num)

    button.clicked.connect(switch_axis)
    refresh(channels_num)
    data_min, data_max = float(np.nanmin(data)), float(np.nanmax(data))
    data_range = data_max - data_min
    cmap = pg.colormap.get(cmap, source="matplotlib")
    bar = pg.ColorBarItem(colorMap=cmap, values=(vmin, vmax), label=cbar_label,
                          interactive=True, rounding=0.0001*data_range)
    bar.setImageItem(image, insert_in=plot)
    h_layout.addWidget(scale_button := get_colorscale_button())
    scale_button.clicked.connect(lambda: reset_scale(bar, vmin, vmax))
    h_layout.addWidget(get_vminmax_button(bar))

    proxy_container = QtWidgets.QGraphicsProxyWidget()
    proxy_container.setWidget(container)
    win.addItem(proxy_container, row=2, col=0)
    pg.exec()

"""Helper Functions"""

def get_layout(size: tuple = (1200, 600), win_title: str = None,
               x_is_time: bool = False):

    """Returns common basic layout of plots."""

    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication(sys.argv)

    win = pg.GraphicsLayoutWidget(show=True)
    win.setWindowTitle(f"FoBench: {win_title}")
    win.setWindowIcon(QtGui.QIcon(str(Path(__file__).resolve().parent/"logo.png")))
    win.setBackground("w")
    win.resize(*size)
    plot = win.addPlot(title=win_title)
    plot.setTitle(win_title, size="18pt", color="k")

    if x_is_time:
        plot.setAxisItems({"bottom": pg.DateAxisItem(utcOffset=1)})
    for axis in (y_axis := plot.getAxis("left"), x_axis := plot.getAxis("bottom")):
        axis.setPen(pg.mkPen("k", width=2))
        axis.setTextPen(pg.mkPen("k"))
        axis.setStyle(tickFont=pg.Qt.QtGui.QFont("Arial", 14))
    return win, app, plot, y_axis, x_axis

def tracker_factory(plot: pg.PlotItem, label: pg.LabelItem, label_text: str,
                    dt: float = None, dx: float = None, dy: float = None,
                    item: pg.ImageItem = None):

    """Factory that returns a function tracking the mouse position and displays it in data units."""

    x_step = dt or dx

    def mouse_moved(evt):
        pos = evt[0]
        if not plot.sceneBoundingRect().contains(pos):
            label.setText("")
            return
        mp = plot.vb.mapSceneToView(pos)
        x = round(mp.x()/x_step)*x_step
        y = round(mp.y()/dy)*dy if dy is not None else mp.y()
        if dt is not None:
            x = datetime.datetime.utcfromtimestamp(x).strftime("%Y-%m-%d %H:%M:%S.%f")
        if item is not None and item.image is not None:
            mp2 = plot.vb.mapFromViewToItem(item, mp)
            col = round(mp2.x())
            row = round(mp2.y())
            if 0 <= col < item.image.shape[0] and 0 <= row < item.image.shape[1]:
                z = item.image[col, row]
                label.setText(label_text.format(x=x, y=y, z=z))
            else:
                label.setText(label_text.rsplit(" | ", 1)[0].format(x=x, y=y))
        else:
            label.setText(label_text.format(x=x, y=y))

    return mouse_moved

def get_axis_button() -> QtWidgets.QPushButton:
    """Returns a button that can be used for switching axis between optical distance and channel number."""
    button = QtWidgets.QPushButton("Distance")
    button.setToolTip("Switch x Axis between Channel Number and Optical Distance")
    button.setFixedWidth(70)
    return button

def get_colorscale_button() -> QtWidgets.QPushButton:
    """Returns a button that can be used for resetting colorscale."""
    button = QtWidgets.QPushButton("Reset colorscale")
    button.setToolTip("Resets the colorscale to initial values")
    button.setFixedWidth(110)
    return button

def get_vminmax_button(bar: pg.ColorBarItem) -> QtWidgets.QPushButton:
    """Returns a button to open dialog to set vmin/vmax values."""
    button = QtWidgets.QPushButton("Set Range")
    button.setToolTip("Modify vmin and vmax of colorscale")
    button.setFixedWidth(70)

    def open_vminmax_dialog():
        dialog = QtWidgets.QDialog()
        dialog.setWindowTitle("Set Colorbar Range")
        dialog.setFixedWidth(280)
        layout = QtWidgets.QFormLayout(dialog)
        layout.setSpacing(10)

        o_vmin, o_vmax = bar.levels()
        vmin_input = QtWidgets.QLineEdit(f"{o_vmin:.6g}")
        vmax_input = QtWidgets.QLineEdit(f"{o_vmax:.6g}")
        layout.addRow("vmin:", vmin_input)
        layout.addRow("vmax:", vmax_input)

        error_label = QtWidgets.QLabel("")
        error_label.setStyleSheet("color: red; font-size: 11px;")
        layout.addRow(error_label)

        btn_box = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        layout.addRow(btn_box)

        def apply():
            try:
                new_vmin = float(vmin_input.text())
                new_vmax = float(vmax_input.text())
            except ValueError:
                error_label.setText("Please enter valid numbers.")
                return
            if new_vmin >= new_vmax:
                error_label.setText("vmin must be less than vmax.")
                return
            error_label.setText("")
            bar.setLevels((new_vmin, new_vmax))
        vmin_input.textChanged.connect(apply)
        vmax_input.textChanged.connect(apply)

        def on_cancel():
            bar.setLevels((o_vmin, o_vmax))
            dialog.reject()

        btn_box.accepted.connect(dialog.accept)
        btn_box.rejected.connect(on_cancel)
        dialog.exec_()
    button.clicked.connect(open_vminmax_dialog)

    return button

def reset_scale(bar: pg.ColorBarItem, vmin: float, vmax: float) -> None:
    """Function to resest the colorscale to original values."""
    bar.setLevels(values=(vmin, vmax))

def get_bottom_layout():
    """Returns bottom layout that holds label for tracking mouse."""
    container = QtWidgets.QWidget()
    h_layout = QtWidgets.QHBoxLayout(container)
    h_layout.setContentsMargins(5, 5, 5, 5)
    text_label = QtWidgets.QLabel()
    text_label.setStyleSheet("color: black; font-size: 10pt;")
    h_layout.addWidget(text_label)
    return text_label, h_layout, container

def get_colors(n: int, colormap:str = "tab10") -> list[tuple[int, int, int]]:
    """Returns a set of n unique colors from a plt colormap."""
    cmap = plt.get_cmap(colormap, n)
    return [tuple(int(x*255) for x in cmap(i)[:3]) for i in range(n)]

"""Special Plotting Functions"""

def plot_fk(wf_ini: np.ndarray, wf_filt: np.ndarray, wf_fk: np.ndarray,
            mask: np.ndarray, f: np.ndarray, k: np.ndarray,
            dt: float) -> None:

    """Plots the in- and outputs of fk filter: initial and filtered wavefield,
    the fk spectrum and the fk mask.

    Parameters
    ----------
    wf_ini : np.ndarray
        Unfiltered wavefield.
    wf_filt : np.ndarray
        Filtered Wavefield.
    wf_fk : np.ndarray
        fk spectrum.
    mask : np.ndarray
        fk mask.
    f : np.ndarray
        Frequency vector.
    k : np.ndarray
        Wavenumber vector.
    dt : float
        Sampling period of data.

    Returns
    -------
    None

    """

    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication(sys.argv)
    win = pg.GraphicsLayoutWidget(show=True)
    win.setWindowTitle("FoBench: Frequency-Wavenumber Filter")
    win.setWindowIcon(QtGui.QIcon(str(Path(__file__).resolve().parent / "logo.png")))
    win.setBackground("w")
    win.resize(1200, 600)

    titles = ["Initial Wavefield", "Filtered Wavefield", "|fk| Spectrum", "fk Mask"]
    images = [wf_ini, wf_filt, np.abs(wf_fk), mask]
    cmaps = [pg.colormap.get("seismic", source="matplotlib"), pg.colormap.get("seismic", source="matplotlib"),
             pg.colormap.get("magma", source="matplotlib"), pg.colormap.get("gray", source="matplotlib")]

    n_samples, n_channels  = wf_ini.shape
    tmax = n_samples*dt
    fmin, fmax = f[0], f[-1]
    kmin, kmax = k[0], k[-1]

    plots = []
    for i, (title, data, cmap) in enumerate(zip(titles, images, cmaps)):
        row, col = divmod(i, 2)
        is_wavefield = i < 2

        plot = win.addPlot(row=row, col=col)
        plot.setTitle(title, color="k", size="14pt")
        plots.append(plot)
        plot.setCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))
        plot.setMouseEnabled(x=True, y=True)
        plot.getViewBox().setMouseMode(pg.ViewBox.RectMode)
        plot.setAspectLocked(False)
        img_item = pg.ImageItem(image=data)
        plot.addItem(img_item)

        if title == "fk Mask":
            vmin, vmax = 0, 1
        else:
            vmax = np.percentile(abs(data), 95)
            vmin = 0 if title == '|fk| Spectrum' else -vmax
        data_min, data_max = float(np.nanmin(data)), float(np.nanmax(data))
        data_range = data_max - data_min
        bar = pg.ColorBarItem(colorMap=cmap, values=(vmin, vmax), interactive=True,
                              rounding=0.0001 * data_range)
        bar.setImageItem(img_item, insert_in=plot)

        for axis in (plot.getAxis("left"), plot.getAxis("bottom")):
            axis.setPen(pg.mkPen("k", width=1))
            axis.setTextPen(pg.mkPen("k"))
            axis.setStyle(tickFont=QtGui.QFont("Arial", 10))

        if is_wavefield:
            img_item.setRect(QtCore.QRectF(0, 0, tmax, n_channels))
            plot.setLimits(xMin=0, xMax=tmax, yMin=0, yMax=n_channels,
                           minXRange=dt, minYRange=1,
                           maxXRange=tmax, maxYRange=n_channels)
            plot.setRange(xRange=(0, tmax), yRange=(0, n_channels), padding=0)
        else:
            img_item.setRect(QtCore.QRectF(fmin, kmin, fmax-fmin, kmax-kmin))
            plot.setLimits(xMin=fmin, xMax=fmax, yMin=kmin, yMax=kmax,
                           minXRange=abs(fmax-fmin) * 0.01,
                           minYRange=abs(kmax-kmin) * 0.01,
                           maxXRange=abs(fmax-fmin),
                           maxYRange=abs(kmax-kmin))
            plot.setRange(xRange=(fmin, fmax), yRange=(kmin, kmax), padding=0)

    labels = [("Channel", "Time [s]"), ("Channel", "Time [s]"),
              ("Wavenumber [1/m]", "Frequency [Hz]"), ("Wavenumber [1/m]", "Frequency [Hz]")]
    for plot, (left, bottom) in zip(plots, labels):
        plot.setLabel("left", left)
        plot.setLabel("bottom", bottom)

    plots[1].setXLink(plots[0])
    plots[1].setYLink(plots[0])
    plots[3].setXLink(plots[2])
    plots[3].setYLink(plots[2])

    pg.exec()