"""Functions for plotting on ``Interrogator`` class level.

:Authors:
    - Sergio Diaz-Meza
    - Jonas Pätzel

:Contributors:
    - Christopher Wollin

"""
import datetime as datetime
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle
from matplotlib.dates import num2date

"""Necessary for Technical aspects"""

class PrecisionDateFormatter(ticker.Formatter):
    """Extends the :class:`~matplotlib.ticker.Formatter` class to allow for millisecond
    precision when formatting a tick (in days since the epoch) with a
    `~datetime.datetime.strftime` format string. Adapted from StackOverflow.
    """

    def __init__(self, fmt:str , precision: int = 3, tz=None):
        """
        Parameters
        ----------
        fmt : str
            `~datetime.datetime.strftime` format string.
        precision : int, optional
            Number of digits for millisecond precision, by default 3.
        tz : timezone, optional
            Timezone for the dates, default is UTC.

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

class AdaptivePrecisionDateFormatter(ticker.Formatter):
    """Extends the :class:`~matplotlib.ticker.Formatter` class to allow for adaptive
    formatting based on the date range, adding millisecond precision if needed.
    """

    def __init__(self, dates, precision=3, tz=None):
        """
        Parameters
        ----------
        dates : array-like
            A list or array of dates to determine the range.
        precision : int, optional
            Number of digits for millisecond precision, by default 3.
        tz : timezone, optional
            Timezone for the dates, default is UTC.

        """
        self.num2date = num2date
        self.dates = dates
        self.precision = precision
        self.tz = tz if tz is not None else datetime.timezone.utc
        self.axis_fmt, self.title_fmt = self._determine_formats()

    def _determine_formats(self):
        """Determine x-axis and title formats based on date range."""
        date_min, date_max = min(self.dates), max(self.dates)
        date_span = date_max - date_min

        if date_span > 30:  # more than a month
            axis_fmt = "%Y-%m-%d"
            title_fmt = ""
        elif date_span > 1:  # Up to a month but more than a day
            axis_fmt = "%d %H:%M"
            title_fmt = "%Y-%m"
        elif date_span <= 1:  # Less than a day
            axis_fmt = "%H:%M:%S.{ms}"
            title_fmt = "%Y-%m-%d"
        else:
            axis_fmt = "%H:%M:%S"
            title_fmt = "%Y-%m-%d %H:%M"

        return axis_fmt, title_fmt

    def __call__(self, x, pos=0):
        # Convert x-axis position to a datetime object
        dt = self.num2date(x, self.tz)

        # Apply milliseconds only if the format requires it
        if "{ms}" in self.axis_fmt:
            ms = dt.strftime("%f")[:self.precision]
            return dt.strftime(self.axis_fmt).format(ms=ms)
        else:
            return dt.strftime(self.axis_fmt)

    def get_formats(self):
        """Return the formats for the x-axis and title."""
        title = self.num2date(self.dates[0], self.tz)

        return self.axis_fmt, title.strftime(self.title_fmt)

"""Interrogator class plots"""

def plot_data_coverage(dataset_infos: str, min_max_t: tuple = None):
    """Plots the coverage of Datasets from a ``Interrogator`` instance.

    Parameters
    ----------
    dataset_infos : str
        Dataframe indicating paths and essential metadata from each file.
    min_max_t : tuple, optional
        Tuple containing the earliest and latest usage time of the Interrogator.

    Returns
    -------
    None

    """

    fig, ax = plt.subplots(1,1, figsize=(8,4))

    y_ticks, y_labels = [], []
    low_lim, high_lim = min_max_t[0].matplotlib_date, min_max_t[1].matplotlib_date # convert to matplotlib datetimes

    for i in range(len(dataset_infos)): # get essential info from datasets.

        # ax.plot()
        dataset_info = dataset_infos[i]# single info line of the dataset.
        low, high = i, i+1
        i_time, f_time = dataset_info # for the moment is just start and end time. More to come! :)
        i_time, f_time = i_time.matplotlib_date, f_time.matplotlib_date

        y_ticks.append(i + 0.5)
        y_labels.append(f'Dataset {high}')

        width, height = f_time - i_time, high - low # calculating dimensions of the rectangle pathc.
        ax.add_patch(Rectangle((i_time, low), width, height)) # Draw the rectangle indicating data availability.

    ax.set_xlim(low_lim, high_lim)
    ax.set_ylim(-1, i+2)

    # Introduce custom ticks and labels
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)

    # Initialize the formatter
    formatter = AdaptivePrecisionDateFormatter([low_lim, high_lim])
    axis_fmt, title_fmt = formatter.get_formats()

    ax.xaxis_date()
    ax.xaxis.set_major_formatter(formatter)
    # precision = str(datetime.timedelta(days=(high_lim - low_lim)).total_seconds())[::-1].find('.')
    # ax.xaxis.set_major_formatter( PrecisionDateFormatter('%H:%M:%S.{ms}', precision) )

    ax.set_xlabel('Time', fontsize=15)
    ax.set_ylabel('Dataset', fontsize=15)

    # plt.xlabel('Time', fontsize=15)
    # plt.ylabel('Datasets', fontsize=15)
    plt.title(title_fmt, fontsize=15)

    plt.tight_layout()
    plt.show()