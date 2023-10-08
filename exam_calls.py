import altair as alt
import pandas as pd
import numpy as np
import scipy.stats as st
from glob import glob
from datetime import date

##### make calls #####
def call_MRD_threshold(input_df, estimated_value_name, output_column="mrd_call", threshold_column="threshold"):
    cancer_type_list = input_df.cancer_type.unique().tolist()
    output_df_list = []
    for CT in cancer_type_list:
        temp_cancer_df = input_df.query(f"cancer_type=='{CT}'").reset_index(drop=True)
        temp_threshold = temp_cancer_df.query(
            "control_group == 'neg'"
        )[estimated_value_name].quantile(0.95, interpolation="higher")
        # drop existing output columns
        drop_columns = [
            col_name for col_name in [output_column, threshold_column] if 
            col_name in temp_cancer_df.columns
        ]
        temp_cancer_df = temp_cancer_df.drop(columns=drop_columns)
        # make calls based on threshold
        temp_cancer_df = temp_cancer_df.assign(
            call_column = temp_cancer_df[estimated_value_name]>temp_threshold,
            thre_column = temp_threshold
        ).rename(columns={"call_column":output_column, "thre_column":threshold_column})
        output_df_list.append(temp_cancer_df)
    return pd.concat(output_df_list, ignore_index=True)


def call_MRD_threshold_with_random_normals(merged_theta_df, genomic_col="ti_genomic_pp"):
    output_list = []
    for random_group_id in merged_theta_df.normal_selection.unique():
        temp_merged_theta_df = merged_theta_df.query(
            f"normal_selection.isna()|normal_selection=='{random_group_id}'"
        )
        temp_combined_call = call_MRD_threshold(temp_merged_theta_df, genomic_col, "genomic_call", "genomic_thre")
        temp_combined_call = call_MRD_threshold(temp_combined_call, "theta_agg_exp", "mbd_call", "mbd_thre")
        temp_combined_call = call_MRD_threshold(temp_combined_call, "combined_theta", "mrd_call", "combined_thre")
        output_list.append(temp_combined_call.assign(normal_selection=random_group_id))
    result_df = pd.concat(output_list, ignore_index=True)
    return result_df

def average_ss_tables(input_df, call_col, subtype_cols):
    output_list = []
    for i in range(1, 101):
        random_group_id = f"group_{i}"
        temp_input_df = input_df.query(
            f"normal_selection=='{random_group_id}'"
        )
        output_list.append(
            call_to_ss_table_anysubtype(temp_input_df, call_col, subtype_cols).assign(
                normal_selection = random_group_id
            )
        )
    concat_df = pd.concat(output_list, ignore_index=True)
#     return concat_df
    result_mean_df = concat_df.groupby(["cancer_type"]+subtype_cols).mean()[[
        "TP", "FN", "FP", "TN", "pos_size", "sensitivity", "specificity"
    ]]
    result_sd_df = concat_df.groupby(["cancer_type"]+subtype_cols).std()[["sensitivity", "specificity"]]
    result_df = result_mean_df.merge(result_sd_df, suffixes=["_mean", "_std"], left_index=True, right_index=True)
    return result_df.reset_index()

##### make plots #####

def plot_estimated_values(input_df, y_column, color_column, threshold_column, cancer_type=None):
    if cancer_type is None:
        temp_plot_data = input_df
    else:
        temp_plot_data = input_df.query(f"cancer_type=='{cancer_type}'")
    base = alt.Chart()
    chart = alt.layer(
        base.mark_point().encode(
            x=alt.X("run_sample_id_plasma", sort='y'),
            y=alt.Y(f"{y_column}:Q"),
            color=f"{color_column}"
        ),
        base.mark_rule(color="red").encode(y=f"{threshold_column}:Q"),
        data = temp_plot_data
    ).properties(width=alt.Step(4)).facet(column="stage_group:N").resolve_scale(x="independent")
    return chart

def plot_estimated_values_log(input_df, y_column, color_column, threshold_column, cancer_type=None):
    if cancer_type is None:
        temp_plot_data = input_df
    else:
        temp_plot_data = input_df.query(f"cancer_type=='{cancer_type}'")
    base = alt.Chart()
    chart = alt.layer(
        base.mark_point().encode(
            x=alt.X("run_sample_id_plasma", sort='y'),
            y=alt.Y(f"{y_column}:Q", scale=alt.Scale(type="log")),
            color=f"{color_column}"
        ),
        base.mark_rule(color="red").encode(y=f"{threshold_column}:Q"),
        data = temp_plot_data
    ).properties(width=alt.Step(4)).facet(column="stage_group:N").resolve_scale(x="independent")
    return chart

def plot_compare_sstables(**kargs):
    methods = kargs.keys()
    for m in methods:
        kargs[m] = kargs[m].assign(method=m.replace('_', ' '))
    plot_data = pd.concat(kargs.values(), ignore_index=True)
    plot_data = plot_data.assign(
        sensitivity = (plot_data.sensitivity*100).round(1),
        specificity = (plot_data.specificity*100).round(1)
    ).replace(
        ["stage_i", "stage_ii", "stage_iii"],
        ["Stage I", "Stage II", "Stage III"]
    )
    result_chart = alt.hconcat()
#     for ct in plot_data.cancer_type:
    for ct in plot_data.cancer_type.unique():
        if plot_data.query(f"cancer_type=='{ct}'").empty:
            continue
        base = alt.Chart().encode(
            x = alt.X("method:N", title=None, sort=None)
            ,y = alt.Y("sensitivity:Q", title="Sensitivity (front) / specificity (background) %")
        )
        temp_chart = alt.layer(
            base.mark_bar(opacity=0.15, color="gray").encode(y="specificity:Q"),
            base.mark_bar().encode(
                color = alt.Color(
                    "stage_group:N",
                    legend=alt.Legend(title="Stage")
                )
            ),
            base.mark_text(dy=8, color="white", fontSize=10).encode(text="sensitivity"),
            base.mark_text(dy=8, fontSize=10).encode(y="specificity:Q",text="specificity"),
            data = plot_data.query(f"cancer_type=='{ct}'")
        ).properties(width=alt.Step(28), height=400)\
        .facet(column=alt.Column("stage_group", title=None, header=alt.Header(labelFontSize=12)))\
        .resolve_scale(x="independent").properties(title=ct.capitalize())
        result_chart |= temp_chart
    return result_chart

def plot_ti_sample_meta(sample_df):
    plot_data = sample_df.replace(
        ["stage_i", "stage_ii", "stage_iii", "grade_6", np.nan, "breast_er_her2", "breast_tp"],
        ["Stage I", "Stage II", "Stage III", "NA", "NA", "breast_her2", "breast_her2"]
    ).rename(columns = {"stage_group":"Stage group"})
    plot_data.cancer_type = plot_data.cancer_type.str.replace("_", " ").str.capitalize()
    base = alt.Chart().encode(
        x=alt.X("Stage group:N"),
        y=alt.Y("count():Q"),
    )
    chart = alt.layer(
        base.mark_bar().encode(color="cancer_subtype:N"),
        base.mark_text(dy=-10).encode(text="count()"),
        data = plot_data
    ).properties(width = alt.Step(30)).facet(
        "cancer_type:N"
    ).resolve_scale(x="independent")
    return chart


def call_to_ss_table_anysubtype(call_df, call_col_name, subtype_cols):
    ss_table_stage_list = []
    for cancer_type in call_df.cancer_type.unique():
        temp_matched_df = call_df.query(f"cancer_type=='{cancer_type}'").replace(
            ["stage_i", "stage_ii", "stage_iii", "grade_6", np.nan, "breast_er_her2", "breast_tp"],
            ["Stage I", "Stage II", "Stage III", "NA", "NA", "breast_her2", "breast_her2"]
        )
        
        TP = temp_matched_df.query(
            f"control_group == 'pos' & {call_col_name}==True"
        ).groupby(["cancer_type"]+subtype_cols).count()["patient_id"].rename("TP").reset_index()
        FP = temp_matched_df.query(
            f"control_group == 'neg' & {call_col_name}==True"
        ).groupby(["cancer_type"]).count()["patient_id"].rename("FP").reset_index()
        TN = temp_matched_df.query(
            f"control_group == 'neg' & {call_col_name}==False"
        ).groupby(["cancer_type"]).count()["patient_id"].rename("TN").reset_index()
        FN = temp_matched_df.query(
            f"control_group == 'pos' & {call_col_name}==False"
        ).groupby(["cancer_type"]+subtype_cols).count()["patient_id"].rename("FN").reset_index()

        temp_ss = TP.merge(FN, how="outer").merge(
            FP.merge(TN, how="outer"), on="cancer_type", how="outer"
        ).replace(np.nan, 0).reset_index(drop=True)

        temp_ss = temp_ss.assign(
            pos_size=temp_ss.TP+temp_ss.FN,
            specificity=temp_ss.TN/(temp_ss.FP+temp_ss.TN),
            sensitivity=temp_ss.TP/(temp_ss.TP+temp_ss.FN),
            FPR=temp_ss.FP/(temp_ss.FP+temp_ss.TN),
            cancer_type=cancer_type
        )
        ss_table_stage_list.append(temp_ss)
    return pd.concat(ss_table_stage_list, ignore_index=True)

def output_compare_sensitive_table(data_df, genomic_col, mbd_col, combined_col, subtype_cols):
    genomic_table = call_to_ss_table_anysubtype(data_df, genomic_col, subtype_cols)
    mbd_table = call_to_ss_table_anysubtype(data_df, mbd_col, subtype_cols)
    combined_table = call_to_ss_table_anysubtype(data_df, combined_col, subtype_cols)
    
    stage_df = genomic_table[["cancer_type"]+subtype_cols].assign(
        total_n = (genomic_table.TP+genomic_table.FN).astype(int)
    )
    merge_cols = ["cancer_type", "sensitivity"]+subtype_cols
    output_df = stage_df.merge(
        genomic_table[merge_cols].rename(columns={"sensitivity":"sensitivity_genomic"}),
        on = ["cancer_type"]+subtype_cols, how="outer"
    ).merge(
        mbd_table[merge_cols].rename(columns={"sensitivity":"sensitivity_mbd"}),
        on = ["cancer_type"]+subtype_cols, how="outer"
    ).merge(
        combined_table[merge_cols].rename(columns={"sensitivity":"sensitivity_mrd"}),
        on = ["cancer_type"]+subtype_cols, how="outer"
    )
    sensitivity_cols = ["sensitivity_genomic", "sensitivity_mbd", "sensitivity_mrd"]
    output_df.loc[:,sensitivity_cols] = (output_df.loc[:,sensitivity_cols]*100).round({
        col:1 for col in sensitivity_cols
    })
    output_df = output_df.sort_values(["cancer_type"]+subtype_cols, ignore_index=True)
    return output_df

def output_compare_sensitive_mean(data_df, genomic_col, mbd_col, combined_col, subtype_cols):
    genomic_table = average_ss_tables(data_df, genomic_col, subtype_cols)
    mbd_table = average_ss_tables(data_df, mbd_col, subtype_cols)
    combined_table = average_ss_tables(data_df, combined_col, subtype_cols)
    
    stage_df = genomic_table[["cancer_type"]+subtype_cols].assign(
        total_n = (genomic_table.TP+genomic_table.FN).astype(int)
    )
    merge_cols = ["cancer_type", "sensitivity_mean"]+subtype_cols
    output_df = stage_df.merge(
        genomic_table[merge_cols].rename(columns={"sensitivity_mean":"sensitivity_genomic"}),
        on = ["cancer_type"]+subtype_cols, how="outer"
    ).merge(
        mbd_table[merge_cols].rename(columns={"sensitivity_mean":"sensitivity_mbd"}),
        on = ["cancer_type"]+subtype_cols, how="outer"
    ).merge(
        combined_table[merge_cols].rename(columns={"sensitivity_mean":"sensitivity_mrd"}),
        on = ["cancer_type"]+subtype_cols, how="outer"
    )
    sensitivity_cols = ["sensitivity_genomic", "sensitivity_mbd", "sensitivity_mrd"]
    output_df.loc[:,sensitivity_cols] = (output_df.loc[:,sensitivity_cols]*100).round({
        col:1 for col in sensitivity_cols
    })
    output_df = output_df.sort_values(["cancer_type"]+subtype_cols, ignore_index=True)
    return output_df

def output_compare_sensitive_mean_std(data_df, genomic_col, mbd_col, combined_col, subtype_cols):
    genomic_table = average_ss_tables(data_df, genomic_col, subtype_cols)
    mbd_table = average_ss_tables(data_df, mbd_col, subtype_cols)
    combined_table = average_ss_tables(data_df, combined_col, subtype_cols)
    
    stage_df = genomic_table[["cancer_type"]+subtype_cols].assign(
        total_n = (genomic_table.TP+genomic_table.FN).astype(int)
    )
    merge_cols = ["cancer_type", "sensitivity_mean", "sensitivity_std"]+subtype_cols
    output_df = stage_df.merge(
        genomic_table[merge_cols].rename(columns={"sensitivity_mean":"sensitivity_genomic", "sensitivity_std":"std_genomic"}),
        on = ["cancer_type"]+subtype_cols, how="outer"
    ).merge(
        mbd_table[merge_cols].rename(columns={"sensitivity_mean":"sensitivity_mbd", "sensitivity_std":"std_mbd"}),
        on = ["cancer_type"]+subtype_cols, how="outer"
    ).merge(
        combined_table[merge_cols].rename(columns={"sensitivity_mean":"sensitivity_mrd", "sensitivity_std":"std_mrd"}),
        on = ["cancer_type"]+subtype_cols, how="outer"
    )
    sensitivity_cols = ["sensitivity_genomic", "sensitivity_mbd", "sensitivity_mrd", "std_genomic", "std_mbd", "std_mrd"]
    output_df.loc[:,sensitivity_cols] = (output_df.loc[:,sensitivity_cols]*100).round({
        col:1 for col in sensitivity_cols
    })
    output_df = output_df.sort_values(["cancer_type"]+subtype_cols, ignore_index=True)
    return output_df

def plot_compare_sstables_anysubtype(subtype_col, **kargs):
    methods = kargs.keys()
    for m in methods:
        kargs[m] = kargs[m].assign(method=m.replace('_', ' '))
    plot_data = pd.concat(kargs.values(), ignore_index=True)
    plot_data = plot_data.assign(
        sensitivity = (plot_data.sensitivity*100).round(1),
        specificity = (plot_data.specificity*100).round(1)
    ).replace(
        ["stage_i", "stage_ii", "stage_iii"],
        ["Stage I", "Stage II", "Stage III"]
    )
    plot_data.loc[:,subtype_col] = plot_data[subtype_col]+"; N="+plot_data.pos_size.astype(int).astype(str)
    
    result_chart = alt.hconcat()
#     for ct in plot_data.cancer_type:
    for ct in plot_data.cancer_type.unique():
        if plot_data.query(f"cancer_type=='{ct}'").empty:
            continue
        base = alt.Chart().encode(
            x = alt.X("method:N", title=None, sort=None)
            ,y = alt.Y("sensitivity:Q", title="Sensitivity (front) / specificity (background) %")
        )
        temp_chart = alt.layer(
            base.mark_bar(opacity=0.15, color="gray").encode(y="specificity:Q"),
            base.mark_bar().encode(
                color = alt.Color(
                    f"{subtype_col}:N",
                    legend=alt.Legend(title="Stage")
                )
            ),
            base.mark_text(dy=8, color="white", fontSize=10).encode(text="sensitivity"),
            base.mark_text(dy=8, fontSize=10).encode(y="specificity:Q",text="specificity"),
            data = plot_data.query(f"cancer_type=='{ct}'")
        ).properties(width=alt.Step(28), height=400)\
        .facet(column=alt.Column(subtype_col, title=None, header=alt.Header(labelFontSize=12)))\
        .resolve_scale(x="independent").properties(title=ct.capitalize())
        result_chart |= temp_chart
    return result_chart

def plot_compare_sstables_anysubtype_with_errorbar(subtype_col, **kargs):
    methods = kargs.keys()
    for m in methods:
        kargs[m] = kargs[m].assign(method=m.replace('_', ' '))
    plot_data = pd.concat(kargs.values(), ignore_index=True)
    plot_data = plot_data.assign(
        sensitivity = (plot_data.sensitivity_mean*100).round(1),
        specificity = (plot_data.specificity_mean*100).round(1),
        sensitivity_lb = (plot_data.sensitivity_mean-plot_data.sensitivity_std)*100,
        sensitivity_ub = (plot_data.sensitivity_mean+plot_data.sensitivity_std)*100
    ).replace(
        ["stage_i", "stage_ii", "stage_iii"],
        ["Stage I", "Stage II", "Stage III"]
    )
    # count sample size
    plot_data.loc[:,subtype_col] = plot_data[subtype_col]+"; N="+plot_data.pos_size.astype(int).astype(str)
    
    result_chart = alt.hconcat()
#     for ct in plot_data.cancer_type:
    for ct in plot_data.cancer_type.unique():
        if plot_data.query(f"cancer_type=='{ct}'").empty:
            continue
        base = alt.Chart().encode(
            x = alt.X("method:N", title=None, sort=None)
            ,y = alt.Y("sensitivity:Q", title="Sensitivity (front) / specificity (background) %")
        )
        temp_chart = alt.layer(
            base.mark_bar(opacity=0.15, color="gray").encode(y="specificity:Q"),
            base.mark_bar().encode(
                color = alt.Color(
                    f"{subtype_col}:N",
                    legend=alt.Legend(title="Stage")
                )
            ),
            base.mark_errorbar().encode(
                y = alt.Y("sensitivity_lb:Q", title="Sensitivity (front) / specificity (background) %"),
                y2="sensitivity_ub:Q"
            ),
            base.mark_text(dy=8, color="white", fontSize=10).encode(text="sensitivity"),
            base.mark_text(dy=8, fontSize=10).encode(y="specificity:Q",text="specificity"),
            data = plot_data.query(f"cancer_type=='{ct}'")
        ).properties(width=alt.Step(28), height=400)\
        .facet(column=alt.Column(subtype_col, title=None, header=alt.Header(labelFontSize=12)))\
        .resolve_scale(x="independent").properties(title=ct.capitalize())
        result_chart |= temp_chart
    return result_chart


def plot_compare_sstables_anysubtype_with_errorbar_no_size(subtype_col, **kargs):
    methods = kargs.keys()
    for m in methods:
        kargs[m] = kargs[m].assign(method=m.replace('_', ' '))
    plot_data = pd.concat(kargs.values(), ignore_index=True)
    plot_data = plot_data.assign(
        sensitivity = (plot_data.sensitivity_mean*100).round(1),
        specificity = (plot_data.specificity_mean*100).round(1),
        sensitivity_lb = (plot_data.sensitivity_mean-plot_data.sensitivity_std)*100,
        sensitivity_ub = (plot_data.sensitivity_mean+plot_data.sensitivity_std)*100
    ).replace(
        ["stage_i", "stage_ii", "stage_iii"],
        ["Stage I", "Stage II", "Stage III"]
    )
    # count sample size
#     plot_data.loc[:,subtype_col] = plot_data[subtype_col]+"; N="+plot_data.pos_size.astype(int).astype(str)
    
    result_chart = alt.hconcat()
#     for ct in plot_data.cancer_type:
    for ct in plot_data.cancer_type.unique():
        if plot_data.query(f"cancer_type=='{ct}'").empty:
            continue
        base = alt.Chart().encode(
            x = alt.X("method:N", title=None, sort=None)
            ,y = alt.Y("sensitivity:Q", title="Sensitivity (front) / specificity (background) %")
        )
        temp_chart = alt.layer(
            base.mark_bar(opacity=0.15, color="gray").encode(y="specificity:Q"),
            base.mark_bar().encode(
                color = alt.Color(
                    f"{subtype_col}:N",
                    legend=alt.Legend(title="Stage")
                )
            ),
            base.mark_errorbar().encode(
                y = alt.Y("sensitivity_lb:Q", title="Sensitivity (front) / specificity (background) %"),
                y2="sensitivity_ub:Q"
            ),
            base.mark_text(dy=8, color="white", fontSize=10).encode(text="sensitivity"),
            base.mark_text(dy=8, fontSize=10).encode(y="specificity:Q",text="specificity"),
            data = plot_data.query(f"cancer_type=='{ct}'")
        ).properties(width=alt.Step(28), height=400)\
        .facet(column=alt.Column(subtype_col, title=None, header=alt.Header(labelFontSize=12)))\
        .resolve_scale(x="independent").properties(title=ct.capitalize())
        result_chart |= temp_chart
    return result_chart