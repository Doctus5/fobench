#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Thu Feb  6 22:16:43 2025

@author: joni
'''

import sys
sys.path.append('/home/joni/Dokumente/GEO4D/Software/fobench/')

# necessary imports for plotting functions
# import sys
import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets, QtGui
from PyQt5 import QtCore
from pathlib import Path

# imports only for testing
from fobench.fiber import Fiber
import numpy as np
import datetime

# %% PyQtGraph example plots and functionality
# import pyqtgraph.examples
# pyqtgraph.examples.run()

#%%

def plot_timeseries(timestamps, data, y_label='', title=''):
    '''
    Co-authors: Jonas Pätzel
    Description: 
        generate generic time series plot using PyQtGraph, ideal for channel plots
    :Params:
        - timestamps(type:numpy): array containing Unix timestamps of data
        - data(type: numpy): array containing data to plot
        - y_label(type: str): y-axis label
        - title(type: str): title of plot
    :Return:
        - NA
    '''

    app = QtWidgets.QApplication(sys.argv)
    win = pg.GraphicsLayoutWidget(show=True)
    win.setWindowTitle(title)
    win.setWindowIcon(QtGui.QIcon(str(Path.cwd() / 'logo.png'))) 
    win.setBackground('w') # white bg
    
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

    sys.exit(app.exec_())

#%%
def plot_2d_timeseries(timestamps, y_ticks, data, y_label='', title='', color_fraction=1):
    '''
    Co-authors: Jonas Pätzel
    Description: 
        generate generic matrix plot where one dimension represents time, e.g. 
        waterfall visualisation of data or spectrograms,PSDs, RMSA ...
    :Params:
        - timestamps(type:numpy): array containing Unix timestamps of data
        - data(type: numpy 2D): array containing data to plot
        - y_label(type: str): y-axis label
        - title(type: str): title of plot
        - color_fraction(type: int): determines min and max value of the colorbar. 
              fraction of the minimum or maximum of the full data, depending on 
              which has the larger absolute value
    :Return:
        - NA
    '''
    
    app = QtWidgets.QApplication(sys.argv)
    win = pg.GraphicsLayoutWidget(show=True)
    win.setWindowTitle(title)
    win.setWindowIcon(QtGui.QIcon(str(Path.cwd() / 'logo.png'))) 
    win.setBackground('w')
    
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
    plot.setAxisItems({'bottom': pg.DateAxisItem(utcOffset=1)}) # when is utcoffset necessary??
    x_axis = plot.getAxis('bottom')
    x_axis.setStyle(tickFont=pg.Qt.QtGui.QFont('Arial', 14))
    x_axis.setTextPen(pg.mkPen('k'))
    x_axis.setPen(pg.mkPen('k', width=2))
    date = datetime.datetime.fromtimestamp(timestamps[0]).strftime('%d.%m.%Y')
    x_axis.setLabel(date, **{'color': 'k', 'font-size': '14pt'})
    plot.setXRange(min(timestamps), max(timestamps), padding=0)

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
            label.setText(f'Time: {x_datetime} | {y_label}: {y_val:.2f}', color='k')
    proxy = pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60, slot=mouse_moved)

    # colorbar and histogram
    cmap = pg.colormap.get('seismic', source='matplotlib')
    histo = pg.HistogramLUTItem()
    histo.setImageItem(image)
    histo.gradient.setColorMap(cmap)
    histo.fillHistogram(fill=True, color='k')
    ticks = histo.gradient.listTicks()
    histo.gradient.removeTick(ticks[1][0])
    histo.gradient.removeTick(ticks[3][0])
    
    level_val = max(abs(data.min()), data.max())
    
    histo.setLevels(color_fraction * -level_val, color_fraction * level_val)
    win.addItem(histo) 

    sys.exit(app.exec_())
#%% TESTING

# das = Fiber('/home/joni/Schreibtisch/Dalvik/data/2024_08_30_02h55m29s_HDAS_2DRawData_Strain.h5', 'aragon')
# das.append_distances()
# das.differentiate()
# ch = 250
# data = das.data[:, ch]

# # plot_timeseries(das.times(time_type='unix'), data, y_label=das.units, title=f'Channel {ch} at Distance {das.distances[ch]} m')

# plot_2d_timeseries(das.times(time_type='unix'), das.distances, das.data, y_label='optical distance [m]', title='data plot', color_fraction=0.1)

