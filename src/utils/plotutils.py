from matplotlib.colors import LinearSegmentedColormap
from pandas import DataFrame

import datetime as dt
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns


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
    
def plotEventDistribution(root_name, ref, freq, pos_min, pos_max):
    # Parameters
    doc_w = 800
    doc_h = 200    
    dna_rel_left = 0.05
    dna_rel_right = 0.95
    dna_rel_top = 0.3
    dna_rel_bottom = 0.9
    
    dna_seq_show = False
    dna_seq_font_size = 8
    dna_stroke_colour = "black"
    dna_stroke_width = 2

    end_label_gap = 0.01
    end_label_font_size = 20

    grid_inc = 10
    grid_stroke_colour = "lightgray"
    grid_stroke_width = 1
    
    grid_label_inc = 10
    grid_label_colour = "lightgray"
    grid_label_font_size = 12

    event_c1 = np.array(mpl.colors.to_rgb("Blue"))
    event_c2 = np.array(mpl.colors.to_rgb("Red"))
    event_max_stroke_width = 2

    if dna_seq_show:
        dna_stroke_width = dna_seq_font_size

    # Checking pos_min and pos_max are different (to prevent divide by zero errors)
    if pos_min == pos_max:
        print("WARNING: Min and max sequence positions must be different")
        return

    # Create document
    datetime_str = dt.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    outname = root_name+"_individual_" + datetime_str + ".svg"
    dwg = svg.Drawing(outname, size = ("%spx" % doc_w, "%spx" % doc_h))

    # Calculating key DNA coordinates
    dna_x1 = doc_w*dna_rel_left
    dna_x2 = doc_w*dna_rel_right
    dna_y1 = doc_h*dna_rel_top
    dna_y2 = doc_h*dna_rel_bottom

    # Adding grid lines
    grid_yc = (dna_y1+dna_y2)/2
    dwg.add(svg.shapes.Line((dna_x1, grid_yc), (dna_x2, grid_yc), stroke=grid_stroke_colour, stroke_width=grid_stroke_width))
    grid_min = grid_inc*math.ceil(pos_min/grid_inc)
    grid_max = grid_inc*math.floor(pos_max/grid_inc)
    if grid_max%grid_inc == 0:
        grid_max += 1
    for grid_pos in range(grid_min, grid_max, grid_inc):
        grid_x = dna_x1 + (dna_x2-dna_x1)*((grid_pos-pos_min)/(pos_max-pos_min))
        grid_t_y = dna_y1+dna_stroke_width/2
        grid_b_y = dna_y2-dna_stroke_width/2
        dwg.add(svg.shapes.Line((grid_x, grid_t_y), (grid_x, grid_b_y), stroke=grid_stroke_colour, stroke_width=grid_stroke_width))

    # Adding grid labels
    grid_label_min = grid_label_inc*math.ceil(pos_min/grid_label_inc)
    grid_label_max = grid_label_inc*math.floor(pos_max/grid_label_inc)
    if grid_label_max%grid_label_inc == 0:
        grid_label_max += 1
    for grid_label_pos in range(grid_label_min, grid_label_max, grid_label_inc):
        grid_label_x = dna_x1 + (dna_x2-dna_x1)*((grid_label_pos-pos_min)/(pos_max-pos_min))
        grid_label_y = dna_y1-dna_stroke_width/2-doc_w*end_label_gap
        rot = "rotate(%i,%i,%i)" % (-90,grid_label_x,grid_label_y)
        dwg.add(svg.text.Text(str(grid_label_pos), insert=(grid_label_x,grid_label_y), transform=rot, style="text-anchor:start; baseline-shift:-50%", font_size=grid_label_font_size, fill=grid_label_colour))

    # Adding cleavage event lines
    total = sum(freq.values())
    max_events = max(freq.values())
    min_events = min(freq.values())
    diff_events = max_events-min_events

    # Preventing divide by zero errors
    if diff_events == 0:
        diff_events = 1

    for (cleavage_site_t, cleavage_site_b) in freq.keys():        
        norm_count = (freq.get((cleavage_site_t, cleavage_site_b))-min_events)/diff_events
        # norm_count = freq.get((cleavage_site_t, cleavage_site_b)/total
        
        # Adding line (adding width 1 to ensure everything is visible)
        event_width = event_max_stroke_width*norm_count+1

        event_t_xc = dna_x1 + (dna_x2-dna_x1)*((cleavage_site_t-pos_min)/(pos_max-pos_min))
        event_t_x1 = event_t_xc-event_width
        event_t_x2 = event_t_xc+event_width
        event_t_y = dna_y1+dna_stroke_width/2

        event_b_xc = dna_x1 + (dna_x2-dna_x1)*((cleavage_site_b-pos_min)/(pos_max-pos_min))
        event_b_x1 = event_b_xc-event_width
        event_b_x2 = event_b_xc+event_width
        event_b_y = dna_y2-dna_stroke_width/2

        col = mpl.colors.to_hex((1-norm_count)*event_c1 + norm_count*event_c2)

        dwg.add(svg.shapes.Polygon(points=[(event_t_x1,event_t_y), (event_t_x2,event_t_y), (event_b_x2,event_b_y), (event_b_x1,event_b_y)], fill=col))

        # Add transparency?

    if dna_seq_show:
        # Draw DNA as text sequence
        dna_pos_increment = (dna_x2-dna_x1)/(pos_max-pos_min)
        ref_rc = ref.reverse_complement()
        for dna_pos in range(pos_min,pos_max+1):
            dna_seq_x = dna_x1+dna_pos_increment*(dna_pos-pos_min)
            dna_seq_y1 = dna_y1+dna_seq_font_size*0.375
            dna_seq_y2 = dna_y2+dna_seq_font_size*0.375
            dwg.add(svg.text.Text(str(ref[dna_pos]), insert=(dna_seq_x, dna_seq_y1), style="text-anchor:middle", font_size=dna_seq_font_size))
            dwg.add(svg.text.Text(str(ref_rc[dna_pos]), insert=(dna_seq_x, dna_seq_y2), style="text-anchor:middle", font_size=dna_seq_font_size))

    else:
        # Draw DNA as lines
        dwg.add(svg.shapes.Line((dna_x1, dna_y1), (dna_x2, dna_y1), stroke=dna_stroke_colour, stroke_width=dna_stroke_width, style="stroke-linecap:square"))
        dwg.add(svg.shapes.Line((dna_x1, dna_y2), (dna_x2, dna_y2), stroke=dna_stroke_colour, stroke_width=dna_stroke_width, style="stroke-linecap:square"))

    # Add DNA end labels
    end_label_x1 = dna_x1-doc_w*end_label_gap
    end_label_x2 = dna_x2+doc_w*end_label_gap
    end_label_y1 = dna_y1+end_label_font_size*0.375
    end_label_y2 = dna_y2+end_label_font_size*0.375
    dwg.add(svg.text.Text("5'", insert=(end_label_x1, end_label_y1), style="text-anchor:end", font_size=end_label_font_size))
    dwg.add(svg.text.Text("3'", insert=(end_label_x2, end_label_y1), style="text-anchor:start", font_size=end_label_font_size))
    dwg.add(svg.text.Text("3'", insert=(end_label_x1, end_label_y2), style="text-anchor:end", font_size=end_label_font_size))
    dwg.add(svg.text.Text("5'", insert=(end_label_x2, end_label_y2), style="text-anchor:start", font_size=end_label_font_size))

    dwg.save()
