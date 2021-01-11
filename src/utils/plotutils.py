from matplotlib.colors import LinearSegmentedColormap
from pandas import DataFrame

import matplotlib as mpl
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
    
def plotEventDistribution(freq):
    # Constants:
    doc_w = 800
    doc_h = 200    
    dna_rel_left = 0.1
    dna_rel_right = 0.9
    dna_rel_top = 0.2
    dna_rel_bottom = 0.8

    pos_min = 1
    pos_max = 2686 # For pUC19, but should be read automatically later

    dna_stroke_colour = "black"
    dna_stroke_width = 2

    end_label_gap = 0.01

    font_size = 20

    max_event_width = 2
    event_c1 = np.array(mpl.colors.to_rgb("Purple"))
    event_c2 = np.array(mpl.colors.to_rgb("Yellow"))

    # Create document
    dwg = svg.Drawing(filename = "test.svg", size = ("%spx" % doc_w, "%spx" % doc_h))

    # Calculating key DNA coordinates
    dna_x1 = doc_w*dna_rel_left
    dna_x2 = doc_w*dna_rel_right
    dna_y1 = doc_h*dna_rel_top-dna_stroke_width/2
    dna_y2 = doc_h*dna_rel_bottom-dna_stroke_width/2

    # Adding cleavage event lines
    total = sum(freq.values())
    max_events = max(freq.values())
    min_events = min(freq.values())
    for (cleavage_site_t, cleavage_site_b) in freq.keys():        
        norm_count = (freq.get((cleavage_site_t, cleavage_site_b))-min_events)/(max_events - min_events)
        # norm_count = freq.get((cleavage_site_t, cleavage_site_b)/total
        
        # Adding line (adding width 1 to ensure everything is visible)
        event_width = max_event_width*norm_count+1

        event_t_xc = dna_x1 + (dna_x2-dna_x1)*((cleavage_site_t-pos_min)/(pos_max-pos_min))
        event_t_x1 = event_t_xc-event_width
        event_t_x2 = event_t_xc+event_width

        event_b_xc = dna_x1 + (dna_x2-dna_x1)*((cleavage_site_b-pos_min)/(pos_max-pos_min))
        event_b_x1 = event_b_xc-event_width
        event_b_x2 = event_b_xc+event_width

        col = mpl.colors.to_hex((1-norm_count)*event_c1 + norm_count*event_c2)

        dwg.add(svg.shapes.Polygon(points=[(event_t_x1,dna_y1), (event_t_x2,dna_y1), (event_b_x2,dna_y2), (event_b_x1,dna_y2)], fill=col))

        # Add transparency?

    # Draw DNA
    dwg.add(svg.shapes.Line((dna_x1, dna_y1), (dna_x2, dna_y1), stroke=dna_stroke_colour, stroke_width=dna_stroke_width))
    dwg.add(svg.shapes.Line((dna_x1, dna_y2), (dna_x2, dna_y2), stroke=dna_stroke_colour, stroke_width=dna_stroke_width))

    # Add DNA end labels
    dwg.add(svg.text.Text("5'", insert=(dna_x1-doc_w*end_label_gap, dna_y1+(dna_stroke_width/2)), style="text-anchor:end; baseline-shift:-50%", font_size=font_size))
    dwg.add(svg.text.Text("3'", insert=(dna_x2+doc_w*end_label_gap, dna_y1+(dna_stroke_width/2)), style="text-anchor:start; baseline-shift:-50%", font_size=font_size))
    dwg.add(svg.text.Text("3'", insert=(dna_x1-doc_w*end_label_gap, dna_y2+(dna_stroke_width/2)), style="text-anchor:end; baseline-shift:-50%", font_size=font_size))
    dwg.add(svg.text.Text("5'", insert=(dna_x2+doc_w*end_label_gap, dna_y2+(dna_stroke_width/2)), style="text-anchor:start; baseline-shift:-50%", font_size=font_size))

    dwg.save()