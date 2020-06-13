from matplotlib.colors import LinearSegmentedColormap
from pandas import DataFrame

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def plotFrequency1D(freq, freq_5p, freq_3p):
    # Creating a Pandas DataFrames for the input dictionaries
    df = DataFrame(data=list(freq.items()), columns=["Sequence", "Events"])
    df_5p = DataFrame(data=list(freq_5p.items()), columns=["Sequence", "Events"])
    df_3p = DataFrame(data=list(freq_3p.items()), columns=["Sequence", "Events"])
    
    df = df.sort_values(by="Events", ascending=False)
    df_5p = df_5p.sort_values(by="Events", ascending=False)
    df_3p = df_3p.sort_values(by="Events", ascending=False)

    fig, axs = plt.subplots(ncols=2, nrows=2)
    gs = axs[0, 0].get_gridspec()
    axs[0, 0].remove()
    axs[1, 0].remove()
    axbig = fig.add_subplot(gs[:, 0])

    # Creating a bar plot of site frequencies
    sns.barplot(x="Events", y="Sequence", data=df, color="slategrey", ax=axbig)

    # Getting colourmap for pie charts
    # cmap = LinearSegmentedColormap.from_list('ne', ['slategrey', 'white'], N=6)
      
    colours = [(0.4392156862745098, 0.5019607843137255, 0.5647058823529412, 1.0),
    (0.5513725490196079, 0.6015686274509804, 0.6517647058823529, 1.0),
    (0.6635294117647059, 0.7011764705882353, 0.7388235294117647, 1.0),
    (0.775686274509804, 0.8007843137254902, 0.8258823529411765, 1.0)
    ]

    axs[0, 1].pie(df_5p["Events"], labels=df_5p["Sequence"], colors=colours)
    axs[1, 1].pie(df_3p["Events"], labels=df_3p["Sequence"], colors=colours)
    
    # Setting axis labels
    axbig.set_title("Dinucleotide frequency")
    axs[0, 1].set_title("5' site frequency")
    axs[1,1].set_title("3' site frequency")

    plt.show()

def plotFrequency2D(labels, freq):
    cmap = LinearSegmentedColormap.from_list('ne', ['#C6CCD2', 'slategrey'], N=256)
    ax = sns.heatmap(freq, cmap=cmap, linewidths=1)
    ax.set_xlabel("Bottom strand dinucleotide")
    ax.set_ylabel("Top strand dinucleotide")
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels,rotation=0) 

    plt.show()
    