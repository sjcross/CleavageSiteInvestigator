from matplotlib.colors import LinearSegmentedColormap
from pandas import DataFrame

import math
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import svgwrite as svg

def plotFrequency1D(freq, freq_5p, freq_3p, show_percentages=True):
    # Creating a Pandas DataFrames for the input dictionaries
    df = DataFrame(data=list(freq.items()), columns=["Sequence", "Events"])
    df_5p = DataFrame(data=list(freq_5p.items()), columns=["Sequence", "Events"])
    df_3p = DataFrame(data=list(freq_3p.items()), columns=["Sequence", "Events"])
    
    df = df.sort_values(by="Events", ascending=False)
    df_5p = df_5p.sort_values(by="Events", ascending=False)
    df_3p = df_3p.sort_values(by="Events", ascending=False)

    # Converting to percentages if necessary
    units = ""
    if show_percentages:
        df.Events = 100*df.Events/df.Events.sum()
        df_5p.Events = 100*df_5p.Events/df_5p.Events.sum()
        df_3p.Events = 100*df_3p.Events/df_3p.Events.sum()
        units = " (%)"

    fig, axs = plt.subplots(ncols=2, nrows=2)
    gs = axs[0, 0].get_gridspec()
    axs[0, 0].remove()
    axs[1, 0].remove()
    axbig = fig.add_subplot(gs[:, 0])

    # Creating a bar plot of site frequencies
    sns.barplot(x="Events", y="Sequence", data=df, color="slategrey", ax=axbig)
    axbig.set(xlabel="Events%s" % units)

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

def plotFrequency2D(labels, freq, show_percentages=True):
    units = ""
    if show_percentages:
        units = " %"
        tot = sum(sum(freq))
        for i in range(len(freq)):
            freq[i] = [round(100*val/tot) for val in freq[i]]

    cmap = LinearSegmentedColormap.from_list('ne', ['#C6CCD2', 'slategrey'], N=256)
    ax = sns.heatmap(freq, cmap=cmap, linewidths=1, cbar_kws={'label': 'Events%s' % units})
    ax.set_xlabel("Bottom strand dinucleotide")
    ax.set_ylabel("Top strand dinucleotide")
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels,rotation=0) 
    
    plt.show()
    
# def plotEventDistribution():
#     # Create document
#     svg_doc = svg.Drawing(filename = "test.svg", size = ("800px", "600px"))

#     svg_doc.add(svg_doc.rect(insert = (0, 0), size = ("200px", "100px"), stroke_width = "1", stroke = "black", fill = "rgb(255,255,0)"))

#     svg_doc.add(svg_doc.text("Hello World", insert = (210, 110)))

#     print(svg_doc.tostring())

#     svg_doc.save()