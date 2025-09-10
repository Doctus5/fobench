-support for Sintela h5 format

-minor changes here and there for readability, clarity e.g. unnecessary if statements replaces by elif

-removed unnecessary commented lines

-removed all unused imports

-interactive_plot method removed as functionality is replaced by .explore method

-Issue #12: removed 'silixa' default for company, added error if user provides no manufacturer

-as distances are directly added in __init__, removed unnecessary check in restrict_channels

-restrict_channels takes two keywords instead of a tuple and sorts them -> das.restrict_channels(ch0, chf)

-Issue #14: fixed as suggested

-added progress bars to more time consuming computations, e.g. coherence and spectral whitening, this requires tqdm to be installed

-much improved logic and readability of __data__ function

-Issue #13 solved:

    -removed attribute 'dataset' completely, it holds only the reference to h5 datasets and prevents .copy() method to work

    -Fiber.base was initialized but then immediatly deleted, instead of holding the reference to the h5 file its now called Fiber.basefile and holds the initial file name as a string


-detrend_signal now takes a single single signal OR a matrix which is detrended channel-wise

-Fiber.corrected is set to True after instrument correction

-envelope was defined twice, now only single definition in signals module

-demeaning got its own function that is called from Fiber.demean and in preprocessing functionality to be more consistent, it takes 1D and 2D data and has "axis" functionality

-Fiber.filter now uses the previoulsy unused function filt_preprocess

-functionality for computaion and plotting of spatial coherence matrix added

-taper_window is replaced by get_tukey_window, now only calling scipy.signal.windows.tukey, this is also called in filt_preprocess
and Fiber.taper which before redundantly implemented tapering again.

-print statement when reading data replaced by progress bar for consistency

###################################################################################################################################

TODO:

-new folder structure, sorts functionality in different categories, longer methods moved outside Fiber class and replaced by 
simple function calls.

-all imports adapted to new folder structure

-added methods to call both plots from the fiber class

-added progress bars to more time consuming tasks

-integrate Johannes fk plot

-taper call taper method Fiber.taper + filter.taper

-filter error message, high = 1.0??

-replace make_plot logic

MAYBE:
-keep DTS references and logic around or unlikely to be implemented?
-reading of terra15 files, special treatment in __data__ really necessary?
-Fiber.data and Fiber.attributes['data'] hold the same reference to the numpy data array, necessary??
