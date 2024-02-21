#import sys
#sys.path.append('/projects/dasatsea/ETNA2019/codes/program_reader/')
import matplotlib.pyplot as plt

########################### BENCH TESTING FOR DAS CODE: FOBENCH #################################

# Hi again!

# He I will explain the basic of the code and basic functions or how to call them so you play with it. 
# I trust you can manage your way through the code. There you can find detailed comments that works as
# documentation and explains what does every method do. I know you have also your own tools in MatLab.
# I hope you can maybe play with both and compare them so you can spot errors, imporvements, etc.

# I highly value any input. You can write me the suggestions here, if you want by typing "# pjousset: " 
# and right next your comments. Remember that the "#" symbol in Python is comments, equivalent to "%" in MatLab.

# I will leave a copy of the code at the side of here, called "fiber.py", which you can enter and even comment
# also. I also encourage this so you can see in more detail the description of each method, inputs, and processing.
# But keep in mind is not the original one, so the computing is actually done with another file elsewhere.
# This is just to have control on the original code, if something happens ;) 
# "In one's disorder one finds order".

# Just to inform, this folder is for you to play around is this one ("philippe_playground"). This folder is a safe
# place for you to test the code, and any changes you produce in the files here will not change the main code or 
# its behavior.

# I will gonna illustrate how is the structure of the fobench code (fiber.py) and from which files depends directly.


#                     fiber.py
#                        |
#                   pyFunctions
#               _________|_________________________
#               |        |          |              |
#           tools.py  filters.py  plotting.py    read_data.py

          
# So basically, to have in mind, the Fiber class sometimes calls functions inside other files to make processes. 
# By doing this, the code doesn't get saturated, and things are more organized. Keep in mind that this is a code
# for basic processing. Complex functions such as the spatial cross-correlation, stacking, SVD, etc. are done
# outside the Fiber Class, but uses the Fiber Class for reading the info and getting relevant configuration values.

# Now let's do now a basic tutorial into the main functions. The other ones you can maybe search in the fobench code.
# Once you agree with the code here, just remember to activate the environment (see README) and execute the code
# by writing "python testing.py" in the Terminal.
# I would suggest that after each line you insert a "exit()" so the code runs until that line and no more. Afterwards
# you can ove the "exit()" down to other lines to let the rest to execute. Another thing you can do is to comment
# which lines you don't want to execute.

print('I am printing something')
#exit()

#################################################################################################################
##                                          CORE METHODS                                                       ##
#################################################################################################################

# Defining the path to the DAS file to read. Currently is able to read 3 types of files (TDMS Silixa, HDF5 Silixa,
# HDF5 Febus, HDF5 ASN, HDF5 QuantX, HDF5 Terra15, Numpy BAM (just testing) ). Therefore you can replace this variable with the one you want. I just
# here put a sugestion.

#file_path = '/projects/dasatsea/ETNA2019/PDN/PDN2019/DASdata/PDN_2019-0.5km-Carina/CarinaP8_Constellation/Etna_UTC_20190912_113158.340.tdms' # Carina Silixa
file_path = '/home/sergioad/Downloads/DAS_ExampleFiles/terra15treble+/20240122_firsttest_strain_rate_UTC-YMD20240124-HMS172633.580_seq_00000000000.hdf5' # Terra15 strain_rate
#file_path = '/home/sergioad/Downloads/DAS_ExampleFiles/terra15treble+/20240122_firsttest_velocity_UTC-YMD20240124-HMS172710.006_seq_00000000000.hdf5' # Terra15 velocity
#file_path = '/home/sergioad/Downloads/DAS_ExampleFiles/OptoDAS/135002.hdf5' # OptoDAS ASN
#file_path = '/home/sergioad/Downloads/DAS_ExampleFiles/OptaSense_QuantX/TestRecord_2024-02-08T115129Z.h5' # OptaSense QuantX

# We load the file with the fobench code. First import the package/code for loading DAS data and then we read.
# The loading part takes several lines because manages different cases for different formats and manufacturers.
# But this is planned to be moved in a separate code that manages this.

from fobench.fiber import Fiber

das = Fiber(file_path, 'terra15') # It is important to set well who is the company or manufacturer of the instrument (letters in lowcaps: quantx, terra15, asn, etc.)!


# The entry variables are the file path, the name of the company. a third one exists which allows loading only
# a certian range, or specific(s) channels (for now only valid for silixa). It is useful for loading only specific prt of the data without
# overloading the PC memory with the rest. However before using this, and the others, it is good to have an
# overwiev of the data. Therefore to plot the data we just use. As in MatLab, you can also interact and do zooms.
# One it is displayed, yoou can zoom-in or out in the image or in the scale to adjust what you want to know.
# The scale is ranged over the maximum value. So if you want to see other features, is better to adjust the color scale by doing zoom on it.

print(das.spatial_interval, das.total_channels)

das.plot()

# Or plot an specific channel.

das.channel_plot(100)

# A quick thing. DAS module is an Object Class, which means that every chanage we introduce to the class will be saved.
# To avoid this there is a function that makes a copy of the Class to preserve the original object.

das2 = das.copy()

# Now that we have the overview, we can decide what to do. We can print the metadata.

das.metadata()

# Print the specific parameters of this recording. I think the variable names explains by themselves.

print(das.sampling_frequency)
print(das.spatial_interval)
print(das.start_time, das.end_time, das.dt)
print(das.total_channels) # number of channels.
print(das.num_points) # number of time samples.
print(das.gauge_length)
print(das.fiber)
print(das.time_length)
print(das.units) # counts? strain-rate?
print(das.data[:10,100]) # let's just display the first 10 values of the channel 100. Better not display the entire thing :)

# We can apply the "instrument correction" which is converting from counts to strain-rate values.
# For other instrument that are not Silixa, ingore this because it will depend on the units on which they came first.

das.instr_correct()
print(das.units) # the units also update from counts to strain-rate.
print(das.data[:10,100])

# We can trim the data for specific times. We can define times in isoformat ('YY-MM-DDTHH:mm:ss.ms'), but for 
# simplicity let us just have 10 seconds after the start of this file.

t0 = das.start_time
tf = t0 + 10
das.trim(t0,tf)
das.plot() # to check that it works.

# We can restrict the channel(s) we want to see or have in the data.

ch0 = 50 # start channel.
chf = 100 # end channel.
das.restrict_channels((ch0,chf))
print(das.total_channels)
das.plot()

# We can also concatenate 2 files to create one single record. This has an advantage, but when you do this recurrently
# you can explode the memory of the PC. The limit goes to the RAM of the PC/Server. Let's try to just concatenate 2
# consecutive files for continuation. You can concatenate with another file, but it will create a gap, obviously, and
# this gap is still information in the memory with "None" items. Since we affected the Fiber Class, it's better that we
# load again the 2 files (change for a path of files you actually have. Be consistent!).

#file_path1 = '/projects/dasatsea/ETNA2019/PDN/PDN2019/DASdata/PDN_2019-0.5km-Carina/CarinaP8_Constellation/Etna_UTC_20190912_113158.340.tdms'
#file_path2 = '/projects/dasatsea/ETNA2019/PDN/PDN2019/DASdata/PDN_2019-0.5km-Carina/CarinaP8_Constellation/Etna_UTC_20190912_113228.340.tdms'
file_path1 = '/home/sergioad/Downloads/DAS_ExampleFiles/terra15treble+/20240122_firsttest_strain_rate_UTC-YMD20240124-HMS172633.580_seq_00000000000.hdf5'
file_path2 = '/home/sergioad/Downloads/DAS_ExampleFiles/terra15treble+/20240122_firsttest_velocity_UTC-YMD20240124-HMS172710.006_seq_00000000000.hdf5'

das1 = Fiber(file_path1, 'terra15').instr_correct() # lets convert the values right from the beginning.
print(das1.start_time, das1.end_time, das1.time_length)
das1.plot() # plot of the das1, without concatenating another file
das2 = Fiber(file_path2, 'terra15').instr_correct()
print(das2.start_time, das2.end_time, das2.time_length)
das2.plot() # the file to concatenate to das1

das1.concatenate(das2, fill_gaps=0) # notice here das2 is being concatenated to das1. Therefore das1 will be modified, not das2.
# Also is important to indicate that gaps must be filled with 0 values in order to avoid computational problems.

print(das1.start_time, das1.end_time, das1.time_length)
das1.plot() # the result of concatenated files.

# Also notice that this method can not work if you decided to concatenate 2 Fiber classes with different channels length.
# I guess this can be also solved by introducing "None" elements where there is no match in term of channel indexes.

#################################################################################################################
##                                          SIGNAL PROCESSING                                                  ##
#################################################################################################################

# I think after seeing this, we can now go to the signal processing routines! So let us load again the data.
 
#file_path = '/projects/dasatsea/ETNA2019/PDN/PDN2019/DASdata/PDN_2019-0.5km-Carina/CarinaP8_Constellation/Etna_UTC_20190912_113158.340.tdms' # Carina Silixa
#file_path = '/home/sergioad/Downloads/DAS_ExampleFiles/terra15treble+/20240122_firsttest_strain_rate_UTC-YMD20240124-HMS172633.580_seq_00000000000.hdf5' # Terra15 strain_rate
file_path = '/home/sergioad/Downloads/DAS_ExampleFiles/terra15treble+/20240122_firsttest_velocity_UTC-YMD20240124-HMS172710.006_seq_00000000000.hdf5' # Terra15 velocity
#file_path = '/home/sergioad/Downloads/DAS_ExampleFiles/OptoDAS/135002.hdf5' # OptoDAS ASN
#file_path = '/home/sergioad/Downloads/DAS_ExampleFiles/OptaSense_QuantX/TestRecord_2024-02-08T115129Z.h5' # OptaSense QuantX

das = Fiber(file_path,'quantx').instr_correct()

das.metadata()

das.plot()
das.channel_plot(100) # visualize the waveform of a random channel

# Now let's try to filter in lowpass, bandpass, and highpass

lf = 50
hf = 200

print('Testing filtering')
das1 = das.copy()
das1.filter('lowpass',lf)
das1.plot()
das1.channel_plot(100)

das1 = das.copy()
das1.filter('bandpass',(lf,hf)) # for bandpass, the variable is a tuple with 2 item: the low and the high frequency.
das1.plot()
das1.channel_plot(100)

das1 = das.copy()
das1.filter('highpass',hf)
das1.plot()
das1.channel_plot(100)

# We can also resample spatially by adding a virtual channel between the actual ones (upsampling) or removing 
# intermediate channels. To resample more, one must execute the code several times.

das1 = das.copy()
print('Total number of channels before upsampling:', das1.total_channels)
das1.spatial_resample('upsampling') # increasing number of channels
print('Total number of channels after upsampling:', das1.total_channels)

das1 = das.copy()
print('Total number of channels before downsampling:', das1.total_channels)
das1.spatial_resample('downsampling') # decreasing number of channels
print('Total number of channels after downsampling:', das1.total_channels)

# We can also decimate the data. For this, we can use Marius Isken filter ('fir-remez') or Javier Quinteros filter
# ('fir235').

das1 = das.copy()
das1.channel_plot(100, show=False)
das1.decimate(200)

t, val = das1.times('matplotlib'), das1.get_data(100) # These are 2 core methods to get a list of times and the values
# of a channel if someone wants to plot or use the values outside the plotting routines of the fobench module. 
plt.plot(t, val, c='red') # in red curve is the decimated data.
plt.show()

# We can also taper, detaper, integrate and differentiate, calculate amplitude root-mean-square, and soon the saturation 
# correction. There is also a function "asStream" that allows to convert the data to a Stream/Trace class from Obspy. 
# This allows then to manipulate tthe signals with Obspy tools, and even save the data in other file formats as mini-seed,
# SAC, etc. Soon a version for Pyrocko will come. 

#################################################################################################################
##                                          PLOTTING FUNCTIONS                                                 ##
#################################################################################################################

# Now lets see some functions for plotting!

# You already saw some of the plots as the general one, and the single channel plot. But we can start here with spectrograms!
# I must confess that I have my doubts in the spectrogram methods. Even if I normalize, seems to not show well.
# We can get a spectrogram over all channels (as the one in Jousset et al., 2022).

das1 = das.copy()
das1.spectrogram()

# Or the spectrogram of a certain channel.

das1.channel_spectrogram(100)

# If you want, you can also plot it with the waveform...

das1.channel_spectrogram(100, trace=True)

# We can also get spectral curves of one or seevral channels.

das1.spectrum([100,110,50]) # 3 channels as example.

# And finally, still in progress, there is the interactive plot, which shows you the waterfall plot and an specific waveform.
# There is also a yellow line marking where the trace is located in the waterplot. You can also change the value of the channel
# and it will update with the new waveform and yellow line. Takes some time depending on the PC.

das1.interactive_plot(100) # go and change the channel in the plot :)

# Aand that is all! Thank you for going through this! I do appreciate the time you could put into this. However this is no the end. 
# This was just the basic overview so you get familiar with the code or how to call functions, run them, and what the methods are
# truly capable of. Therefore I would like to invite you to exploit the functions of the code in any combination, as if you were 
# gonna do routine processing for different kind of purposes.

# You can pass to another file of your own (a blanck Python file), where the idea is that you exploit the code as far as you can!
# Thank you really, and have fun, I guess? :)


