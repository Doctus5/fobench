#Obspy and Pyrocko stuff!
import obspy as ob
from obspy.core import UTCDateTime as UTC
#from pyrocko.plot.automap import Map

#Normal libraries and matplotlib
import os
import json
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.image import imread
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.widgets import TextBox
import matplotlib.patheffects as PathEffects
from matplotlib.dates import DateFormatter, MinuteLocator, num2date
#from matplotlib_scalebar.scalebar import ScaleBar
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime as datetime
#import utm



'''
########################################################################################
Necessary for Technical aspects...
########################################################################################
'''

# #Setting the precision for the labels in the plots. Taken from StackOverflow.
class PrecisionDateFormatter(ticker.Formatter):
    """
    Extend the `matplotlib.ticker.Formatter` class to allow for millisecond
    precision when formatting a tick (in days since the epoch) with a
    `~datetime.datetime.strftime` format string.
    """

    def __init__(self, fmt, precision=3, tz=None):
        """
        Parameters
        ----------
        fmt : str
            `~datetime.datetime.strftime` format string.
        """
        self.num2date = num2date
        self.fmt = fmt
        self.tz = tz if tz is not None else datetime.timezone.utc
        self.precision = precision

    def __call__(self, x, pos=0):
        if x == 0:
            raise ValueError("DateFormatter found a value of x=0, which is "
                             "an illegal date; this usually occurs because "
                             "you have not informed the axis that it is "
                             "plotting dates, e.g., with ax.xaxis_date()")

        dt = self.num2date(x, self.tz)
        ms = dt.strftime("%f")[:self.precision]

        return dt.strftime(self.fmt).format(ms=ms)
    
		
'''
########################################################################################
DAS Class plots below...
########################################################################################
'''

#Plotting the DAS as an image in the Class DAS.		
def gen_DAS_plot(data=None, t=None, channels=None, units_y=None, max_value=None, figsize=None, title=None, cmap='seismic', show=True, file_name=None, where=None, add_data=None, **kwargs):
	
	fig, ax = plt.subplots() if figsize is None else plt.subplots(figsize=figsize)
	fig.autofmt_xdate()
	max_val = max(data.max(), abs(data.min())) if max_value == None else max_value
	cm = ax.imshow(data, vmin=-max_val, vmax=max_val, cmap=cmap, extent=(channels[0],channels[-1]+1,t[-1],t[0]), aspect='auto', interpolation='none', **kwargs)
	plt.colorbar(cm, label=units_y)
	
	#Additional data
	if add_data is not None:
		
		ax.imshow(add_data, cmap='binary', extent=(channels[0],channels[-1]+1,t[-1],t[0]), aspect='auto', alpha=add_data)
	
	ax.yaxis_date()
	precision = str(datetime.timedelta(days=(t[1]-t[0])).total_seconds())[::-1].find('.')
	ax.yaxis.set_major_formatter( PrecisionDateFormatter('%H:%M:%S.{ms}', precision) )
	
	plt.ylabel('Time', fontsize=15)
	plt.xlabel('Channel', fontsize=15)
	plt.title(title, fontsize=20)
	
	ax.tick_params(axis="x", bottom=True, top=True, labelbottom=True, labelleft=True, labeltop=True, rotation=0)
	
	if show == True:
		plt.show()
	if file_name is not None:
		fig.savefig(file_name, transparent=True, bbox_inches='tight', pad_inches = 0)
		
		
#Function for the interactive plot.
def DAS_interactive_plot(data=None, set_channel=None, t=None, channels=None, units_y=None, max_value=None, figsize=None, title=None, cmap='seismic', show=True, file_name=None, where=None, add_data=None, **kwargs):

	#fig, ax = plt.subplots(2,1) if figsize is None else plt.subplots(2,1,figsize=figsize)
	fig, ax = plt.subplots(2,1, sharex=True, gridspec_kw={'hspace': 0.2,'height_ratios':[3,1]}, constrained_layout=True) if figsize is None else plt.subplots(2,1, sharex=True, gridspec_kw={'hspace': 0.2,'height_ratios':[3,1]}, figsize=figsize)
	fig.autofmt_xdate()
	max_val = max(data.max(), abs(data.min())) if max_value == None else max_value
	pos_channel = channels.index(set_channel)
	
	#plot1
	data_signal = data[:,pos_channel]
	data_M = np.fliplr(data).T
	cm = ax[0].imshow(data_M, vmin=-max_val, vmax=max_val, cmap=cmap, extent=(t[0],t[-1],channels[0],channels[-1]+1), aspect='auto', interpolation='none', **kwargs)
	#plt.colorbar(cm, label=units_y)
	ax[0].plot([t[0],t[-1]],[float(set_channel),float(set_channel)],color='yellow')
	ax[0].set_ylabel('Channel')
	ax[0].set_title(title)
	
	#plot2
	ax[1].plot(t, data_signal, c='black', linewidth=0.7, label=str(set_channel).zfill(5))
	ax[1].set_ylim(-max_val-(max_val*0.2), max_val+(max_val*0.2))
	#ax[1].legend(loc=1)
	
	ax[1].xaxis_date()
	precision = str(datetime.timedelta(days=(t[1]-t[0])).total_seconds())[::-1].find('.')
	ax[1].xaxis.set_major_formatter( PrecisionDateFormatter('%H:%M:%S.{ms}', precision) )
	ax[1].set_ylabel(units_y)
	ax[1].ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
	ax[1].set_xlabel('Time')
	
	fig.colorbar(cm, ax=ax[:], label=units_y)
	
	def update(text):
	
		ax[0].lines[0].remove()
		ax[1].lines[0].remove()
		
		pos_channel = channels.index(int(text))
		data_signal = data[:,pos_channel]
		
		ax[0].plot([t[0],t[-1]],[float(text),float(text)],color='yellow')
		ax[1].plot(t, data_signal, c='black', linewidth=0.7, label=str(set_channel).zfill(5))
		
		print('Updated')
	
	#box
	box = plt.axes([0.684,0.357,0.06,0.03])
	textbox = TextBox(box, 'Channel', initial=set_channel, label_pad=0.1)
	textbox.on_submit(update)
	#textbox.on_text_change()
	#textbox.set_val()
	
	if show == True:
		plt.show()
	if file_name is not None:
		fig.savefig(file_name, transparent=True, bbox_inches='tight', pad_inches = 0, dpi=300)


#Function for simple plot of a channel for the DAS class object.
def simple_plot(data, t, channel='', units_y=None, max_value=None, spectrogram=False, show=True, figsize=None, title=None, file_name=None, where=None, **kwargs):

	fig, ax = plt.subplots(1,1, sharex=True, gridspec_kw={'hspace': 0.3}) if figsize is None else plt.subplots(1,1, sharex=True, gridspec_kw={'hspace': 0.3}, figsize=figsize)
	fig.autofmt_xdate()
	
	ax.plot(t, data, c='black', linewidth=0.7, label=str(channel).zfill(5), **kwargs)
	max_val = max(data.max(), abs(data.min())) if max_value == None else max_value
	ax.set_ylim(-max_val-(max_val*0.2), max_val+(max_val*0.2))
	
	ax.xaxis_date()
	precision = str(datetime.timedelta(days=(t[1]-t[0])).total_seconds())[::-1].find('.')
	ax.xaxis.set_major_formatter( PrecisionDateFormatter('%H:%M:%S.{ms}', precision) )
	ax.set_xlim(t[0],t[-1])
	
	ax.legend(loc=1)
	ax.set_ylabel(units_y, fontsize=15)
	ax.set_xlabel('Time', fontsize=15)
	ax.set_title(title, fontsize=20)
	plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))

	if show == True:
		plt.show()
	if file_name is not None:
		fig.savefig(file_name, transparent=True, bbox_inches='tight', pad_inches = 0)
		

#Function for the DAS spectrogram in the class DAS. Channels and the spectrogram matrix must already be computed and passed as an input. This is done in the DAS class under the method spectrogram().
def gen_spectrogram(spec_matrix=None, freqs=None, x=None, max_value=None, units_y=None, figsize=None, cmap='viridis', title=None, show=True, file_name=None, where=None, **kwargs):

	fig, ax = plt.subplots() if figsize is None else plt.subplots(figsize=figsize)
	
	plt.ylabel('Frequency [Hz]', fontsize=15)
	plt.xlabel('Channel', fontsize=15)
	extent = (x[0], x[-1]+1, freqs[0], freqs[-1])

	cm = ax.imshow(spec_matrix, cmap=cmap, vmin=0, vmax=max_value, extent=extent, aspect='auto', interpolation='none', **kwargs)
	plt.colorbar(cm, ax=ax, label=units_y)
	
	ax.set_title(title, fontsize=20)
	ax.tick_params(axis="x", bottom=True, top=True, labelbottom=True, labelleft=True, labeltop=True, rotation=0)
	
	if show == True:
		plt.show()
	if file_name is not None:
		fig.savefig(file_name, transparent=True, bbox_inches='tight', pad_inches = 0)
		

#Functino for plotting spectrogram in time for one channel in class DAS.
def simple_spectrogram(data=None, freq=None, t=None, units_y=None, figsize=None, trace=None, cmap='viridis', title=None, show=None, file_name=None, where=None, freq_lim=None, **kwargs):

	if trace is not None:

		fig, ax = plt.subplots(2,1, sharex=True, gridspec_kw={'hspace': 0.02,'height_ratios':[3,1]}, constrained_layout=True) if figsize is None else plt.subplots(2,1, sharex=True, gridspec_kw={'hspace': 0.2,'height_ratios':[3,1]}, figsize=figsize)
		
		extent = (t[0], t[-1], freq[0], freq[-1])
		cm = ax[0].imshow(data, cmap=cmap, extent=extent, aspect='auto', interpolation='none',**kwargs)
		plt.colorbar(cm, ax=ax[0], label=units_y)
	
		precision = str(datetime.timedelta(days=(t[1]-t[0])).total_seconds())[::-1].find('.')
		ax[0].xaxis.set_major_formatter( PrecisionDateFormatter('%H:%M:%S.{ms}', precision) )
		ax[0].set_xlim(t[0],t[-1])
		ax[0].set_title(title, fontsize=20)
		ax[0].set_ylabel('Frequency [Hz]', fontsize=10)
		
		ax[1].plot(t, trace, c='black', linewidth=0.7)
		
		ax[1].set_ylabel(units_y, fontsize=10)
		ax[1].set_xlabel('Time', fontsize=15)
		
		if freq_lim is not None:
		
			ax[0].set_ylim(freq_lim[0],freq_lim[1])
	
		ax[1].tick_params(axis="x", bottom=True, top=False, labelbottom=True, labelleft=True, labeltop=False, labelrotation=25)
		
		
	else:

		fig, ax = plt.subplots() if figsize is None else plt.subplots(figsize=figsize)
		fig.autofmt_xdate()

		extent = (t[0], t[-1], freq[0], freq[-1])
		cm = ax.imshow(data, cmap=cmap, extent=extent, aspect='auto', interpolation='none', **kwargs)
		plt.colorbar(cm, ax=ax, label=units_y)
	
		precision = str(datetime.timedelta(days=(t[1]-t[0])).total_seconds())[::-1].find('.')
		ax.xaxis.set_major_formatter( PrecisionDateFormatter('%H:%M:%S.{ms}', precision) )
		ax.set_xlim(t[0],t[-1])
	
		plt.title(title, fontsize=20)
		plt.ylabel('Frequency [Hz]', fontsize=15)
		plt.xlabel('Time', fontsize=15)
	
		ax.tick_params(axis="x", bottom=True, top=False, labelbottom=True, labelleft=True, labeltop=False, labelrotation=25)
	
	if show == True:
		plt.show()
	if file_name is not None:
		fig.savefig(file_name, transparent=True, bbox_inches='tight', pad_inches = 0)
		

#Function to plot simple spectrum (1D) for certain channels of DAS class.
def simple_spectrum(spectrums=None, freqs=None, channels=None, y_units=None, legend=True, figsize=None, show=True, 
                    file_name=None, where=None, title=None, **kwargs):

	fig = plt.figure() if figsize is None else plt.figure(figsize=figsize)
	
	plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
	plt.xlabel('Frequency [Hz]', fontsize=15)
	plt.ylabel(y_units, fontsize=15)
	plt.title(title, fontsize=20)
	plt.xlim(freqs[0],freqs[-1])
	#plt.xscale('log')
	plt.ylim(0,spectrums.max() + spectrums.max()/10)
	#plt.xlim(1E-3,100)
	plt.grid()
	
	for i in range(len(spectrums)):
		
		plt.plot(freqs, spectrums[i], label=channels[i], **kwargs)
	
	if legend == True:
	
		plt.legend()
	
	if show == True:
		plt.show()
	if file_name is not None:
		fig.savefig(file_name, transparent=True, bbox_inches='tight', pad_inches = 0)


# Function for plotting channels as record sections
def plot_record_section(signals, t, channels, date):
    
	num_stations = len(channels)

	fig, ax = plt.subplots()
	fig.set_size_inches(10, 8)

	# Calculate the scaling factor based on the maximum absolute value in the signals
	max_val = np.max(np.abs(signals))
	scaling_factor = 1.5 / max_val

	# Plot the signals as traces with further scaled down values
	for i in range(num_stations):
		#y = num_stations - i
		y = i
		ax.plot(t, signals[:, i] * scaling_factor + y, color='black', linewidth=0.7)

	# Add station labels
	ax.set_yticks(np.arange(0, num_stations))
	ax.set_yticklabels([f'Ch {i}' for i in channels])
	ax.invert_yaxis()

	ax.xaxis_date()
	precision = str(datetime.timedelta(days=(t[1]-t[0])).total_seconds())[::-1].find('.')
	ax.xaxis.set_major_formatter( PrecisionDateFormatter('%H:%M:%S.{ms}', precision) )
	ax.set_xlim(t[0],t[-1])

	ax.grid(color='gray', linestyle='--', alpha=0.8)
	ax.set_xlabel('Time (s)', fontsize=15)
	ax.set_ylabel('Station', fontsize=15)
	ax.set_title('Record Section for '+date, fontsize=20)

    # Show the plot
	plt.show()


				

