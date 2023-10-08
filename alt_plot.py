import altair as alt
import pandas as pd
import numpy as np
import math
alt.data_transformers.disable_max_rows()

##### plot functions #####

def alt_rel_plot(data, x, y, color = None, log_scale = False):
    if color is None:
        plot_data = data[[x, y]].round(2).drop_duplicates()
        chart = alt.Chart(data=plot_data).mark_circle().encode(x=f"{x}:Q", y=f"{y}:Q")
    else:
        plot_data = data[[x, y, color]].round(2).drop_duplicates()
        chart = alt.Chart(data=plot_data).mark_circle().encode(
            x=f"{x}:Q", y=f"{y}:Q",
            color = color
        )
    if log_scale:
        plot_data = data[[x, y]].round(3).drop_duplicates()
        plot_data = plot_data.replace(0, np.nan)
        chart = alt.Chart(data=plot_data).mark_circle().encode(
            x=alt.X(f"{x}:Q", scale=alt.Scale(type='log')),
            y=alt.Y(f"{y}:Q", scale=alt.Scale(type='log'))
        )
    return chart

def alt_hist_plot(data, x = None, step = None, bin_n=50, log_scale = False, rename = None):
    # take series from df #
    if x is None:
        plot_data = data
    else:
        plot_data = data[x]

    trait_name=plot_data.name
    
    # generate intervals #
    if step is not None:
        bin_n = int((plot_data.max() - plot_data.min())/step)+1
    
    plot_data = pd.cut(plot_data, bin_n).reset_index().groupby(trait_name).count(
    ).rename(columns={"index":"Count"}).reset_index()
    plot_data = plot_data.assign(mid_point=pd.IntervalIndex(plot_data[trait_name]).mid)[["Count", "mid_point"]]
    plot_data = plot_data.loc[plot_data.Count>0]

    bar_size = int(350/bin_n)
    if rename is None:
        x_name = trait_name
    else:
        x_name = rename
    # make plot #
    if log_scale:
        plot_data_log = plot_data
        plot_data_log.Count.replace(1, 1.5, inplace=True)
        chart = alt.Chart(plot_data_log).mark_bar(size=bar_size).encode(
            x = alt.X("mid_point:Q", axis=alt.Axis(title=x_name)),
            y = alt.Y("Count:Q", scale=alt.Scale(type='log'))
        )
    else:
        chart = alt.Chart(plot_data).mark_bar(size=bar_size).encode(
            x = alt.X("mid_point:Q", axis=alt.Axis(title=x_name)),
            y = "Count"
        )
    return chart

def alt_plot_grid(plot_list, ncol=None, nrow = None):
    total_n = len(plot_list)
    if total_n <1:
        return None
    # get ncol #
    if (ncol is None)and (nrow is None):
        ncol = round(math.sqrt(total_n))
    elif ncol is None:
        ncol = math.ceil(total_n/nrow)
    # add plots to rows #
    plot_rows = []
    temp_plot = plot_list[0]
    for i in range(1, total_n):
        if i%ncol == 0:
            plot_rows.append(temp_plot)
            temp_plot = plot_list[i]
        else:
            temp_plot = temp_plot | plot_list[i]
    plot_rows.append(temp_plot)
    # aggerate plots #
    out_plot = plot_rows[0]
    for i in range(1, len(plot_rows)):
        out_plot = out_plot & plot_rows[i]
    return out_plot

def plot_roc_curve(input_df, truth_col, pos_label, trait_cols, roc_x_scale=[0,1], roc_y_scale=[0,1]):
    plot_data_list=[]
    for t in trait_cols:
        fpr, tpr, thre = sk.metrics.roc_curve(input_df[truth_col], input_df[t], pos_label=pos_label)
        temp_roc_df = pd.DataFrame({"FPR":fpr, "TPR":tpr, "thresholds":thre}).assign(Trait=t)
        plot_data_list.append(temp_roc_df)
    roc_df = pd.concat(plot_data_list, ignore_index=True)
    abline_df = pd.DataFrame({
        "xy_line":[0,1], "Trait":trait_cols[0]
    })
    plot_data = pd.concat([roc_df, abline_df], ignore_index=True)
    chart=alt.layer(
        alt.Chart().mark_line(color="red", opacity=0.3, clip=True).encode(
            x="xy_line:Q",
            y="xy_line:Q"
        ),
        alt.Chart().mark_line(clip=True).encode(
            x=alt.X("FPR:Q", title="FPR", scale=alt.Scale(domain=roc_x_scale)),
            y=alt.Y("TPR:Q", title="TPR", scale=alt.Scale(domain=roc_y_scale)),
            color="Trait:N"
        ),
        data=plot_data
    ).properties(width=400, height=400)
    return chart