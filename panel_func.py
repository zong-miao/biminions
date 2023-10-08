import pandas as pd
import altair as alt

def parse_gtf_file(gtf_file):
    gtf_df = pd.read_csv(
        gtf_file, sep = '\t', compression='gzip', header=None, low_memory = False
        ,names = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
    )
    attributes_df = gtf_df["attributes"].str.split(' *; *', expand = True)
    gene_id = attributes_df.iloc[:,0].str.split(' ').str[-1].str.strip('"')
    transcript_id = attributes_df.iloc[:,1].str.split(' ').str[-1].str.strip('"')
    gtf_df["gene_id"] = gene_id
    gtf_df["transcript_id"] = transcript_id
    
    ## convert 1 base close region to 0 base open ended ##
    gtf_df.loc[:,"start"] = gtf_df["start"]-1
    return gtf_df.drop(columns=["attributes"])

def fetch_bed_file(file_name):
    """
    Fetches the bed annotation from probe bed file,
    seperate the annotaiton part into columns
    """
    df = pd.read_csv(file_name, sep="\t",skiprows=1, header=None,
                     names=["chrom", "start", "end", "ped_probe_id", "score", "annot"])

    # get column names #
    annot_dict = {}
    for annot in list(df["annot"]):
        column_names = {data[0]:[] for data in [x.split('=') for x in annot.split(';')]}
        annot_dict.update(column_names)

    # load data from each row #
    for annot in list(df["annot"]):
        line_dict = {data[0]:data[1] for data in [x.split('=') for x in annot.split(';')]}
        for item in annot_dict:
            if item in line_dict:
                annot_dict[item].append(line_dict[item])
            else:
                annot_dict[item].append("NA")
    annot_df = pd.DataFrame(annot_dict)
    return pd.concat((df.iloc[:,0:5], annot_df), axis=1)

def fetch_gff_file(file_name):
    """
    Fetches the gff annotation from gff file,
    seperate the annotaiton part into columns
    """
    df = pd.read_csv(file_name, sep="\t",skiprows=1, header=None,
                     names=["NAME", "SOURCE", "TYPE", "START", "END", "SCORE", "STRAND", "PHASE", "ANNOT"])

    # get column names #
    annot_dict = {}
    for annot in list(df["ANNOT"]):
        column_names = {data[0]:[] for data in [x.split('=') for x in annot.split(';')]}
        annot_dict.update(column_names)

    # load data from each row #
    for annot in list(df["ANNOT"]):
        line_dict = {data[0]:data[1] for data in [x.split('=') for x in annot.split(';')]}
        for item in annot_dict:
            if item in line_dict:
                annot_dict[item].append(line_dict[item])
            else:
                annot_dict[item].append("NA")
    annot_df = pd.DataFrame(annot_dict)
    return pd.concat((df.iloc[:,0:8], annot_df), axis=1)

def merge_regions(*args):
    input_df = pd.concat(args)
    input_df = input_df.astype({"chrom":str, "start":int, "end":int}).sort_values(["chrom", "start", "end"]).reset_index(drop=True)
    input_df = input_df.drop_duplicates(subset=["chrom", "start", "end"])
    
    # remove the embeded regions: 1-1000 fowllowed by 10:20
    embeded = (
        (input_df["chrom"] == input_df["chrom"].shift()) &
        (input_df["start"] >= input_df["start"].shift()) &
        (input_df["end"] <= input_df["end"].shift())
    )
    while sum(embeded) > 0:
        input_df = input_df.loc[~embeded]
        embeded = (
            (input_df["chrom"] == input_df["chrom"].shift()) &
            (input_df["start"] >= input_df["start"].shift()) &
            (input_df["end"] <= input_df["end"].shift())
        )
        
    # mark overlapped 
    input_df["overlap_group"]=(
        (input_df["start"] > input_df["end"].shift()) |
        (input_df["chrom"] != input_df["chrom"].shift())
    ).cumsum()
    unique_join = lambda x:';'.join(set(x))
    agg_par_dict = {k:"first" for k in input_df.columns if k != "overlap_group"}
    agg_par_dict.update({"start":"min", "end":"max"})
    merged_df = input_df.groupby("overlap_group").agg(agg_par_dict).reset_index(drop=True)
    return merged_df

def find_overlapped_region(input_df, target_region):
    target_region = merge_regions(target_region[["chrom", "start", "end"]])
    summed_df = pd.concat([
        input_df.assign(temp_mark=0),
        target_region.assign(temp_mark=1)
    ]).sort_values(["chrom", "start", "end"], ignore_index=True).astype({"start":int, "end":int})
    
    # get the closest region
    l_target_index = summed_df["temp_mark"].cumsum()-1
    l_target_index.loc[l_target_index<0] = 0
    r_target_index = l_target_index + 1
    r_target_index.loc[r_target_index>=(target_region.shape[0]-1)] = (target_region.shape[0]-1)
    
    left_target_regions = target_region.iloc[l_target_index,:].reset_index(drop=True)
    left_target_regions.columns = ["chrom", "start", "end"]
    right_target_regions = target_region.iloc[r_target_index,:].reset_index(drop=True)
    right_target_regions.columns = ["chrom", "start", "end"]
    
    # mark overlapped regions
    summed_df.loc[:,"temp_mark1"] = (
        (
            (summed_df["chrom"] == left_target_regions["chrom"])
            & (summed_df["start"] < left_target_regions["end"].astype(int))
            & (summed_df["end"] > left_target_regions["start"].astype(int))
        ) | (
            (summed_df["chrom"] == right_target_regions["chrom"])
            & (summed_df["start"] < right_target_regions["end"].astype(int))
            & (summed_df["end"] > right_target_regions["start"].astype(int))
        )
    )
    
    # get overlapped regions
    result = summed_df.loc[summed_df["temp_mark1"] & (summed_df["temp_mark"] == 0)].drop(columns=["temp_mark", "temp_mark1"])
    return result.reset_index(drop=True)

def extract_overlapped_region(input_df, target_region):
    target_region = merge_regions(target_region[["chrom", "start", "end"]])
    target_region = find_overlapped_region(target_region, input_df)
    if target_region.shape[0]<1:
        return target_region
    summed_df = pd.concat([
        input_df.assign(temp_mark=0),
        target_region.assign(temp_mark=1)
    ]).sort_values(["chrom", "start", "end"], ignore_index=True).astype({"start":int, "end":int})
    
    # get the closest region
    l_target_index = summed_df["temp_mark"].cumsum()-1
    l_target_index.loc[l_target_index<0] = 0
    
    def update_position(target_index):
        temp_target_regions = target_region.iloc[target_index,:].reset_index(drop=True).rename(
            columns = {"chrom":"chrom_alt", "start":"start_alt", "end":"end_alt"}
        )
        temp_merged_df = pd.concat([summed_df, temp_target_regions], axis=1)
        temp_merged_df = temp_merged_df.loc[temp_merged_df["temp_mark"] == 0]
        temp_merged_df = temp_merged_df.loc[(
            (temp_merged_df["chrom"] == temp_merged_df["chrom_alt"])
            & (temp_merged_df["start"] < temp_merged_df["end_alt"].astype(int))
            & (temp_merged_df["end"] > temp_merged_df["start_alt"].astype(int))
        )]
        temp_merged_df.loc[:,"start"] = temp_merged_df[["start", "start_alt"]].max(axis=1)
        temp_merged_df.loc[:,"end"] = temp_merged_df[["end", "end_alt"]].min(axis=1)
        return temp_merged_df.astype({"start":int, "end":int})
    
    df_list = []
    # check the nearest left target_index
    df_list.append(update_position(l_target_index))
    
    # check the nearest right target_index until no match
    # this can be optimized
    max_idx = target_region.shape[0]-1
    i=1
    while i > 0 and i <= max_idx:
        r_target_index = l_target_index + i
        r_target_index.loc[r_target_index>=max_idx] = max_idx
        right_overlapped_df = update_position(r_target_index)
        if right_overlapped_df.shape[0] > 0:
            df_list.append(right_overlapped_df)
            i += 1
        else:
            i = 0
    
    # output the overlaped regions
    result = pd.concat(df_list, ignore_index=True)\
    .drop_duplicates(subset=["chrom", "start", "end"]).\
    sort_values(["chrom", "start", "end"]).reset_index(drop=True)

    return result.drop(columns = ["chrom_alt", "start_alt", "end_alt", "temp_mark"])

def complement_regions(target_region, max_posi=1e10):
    target_region = merge_regions(target_region)
    df_list=[]
    for chromosome in target_region["chrom"].unique():
        temp_df = target_region.loc[target_region["chrom"] == chromosome]
        positions = sorted(temp_df["start"].tolist() + temp_df["end"].tolist() + [0, max_posi])
        df_list.append(pd.DataFrame({
            "chrom":chromosome,
            "start":[positions[i] for i in range(0, len(positions), 2)],
            "end":[positions[i] for i in range(1, len(positions), 2)],
        }))
    return pd.concat(df_list, ignore_index=True).astype({"start":int, "end":int})

def extract_non_overlapped_region(input_df, target_region):
    # take complement of the target regions
    target_complement = complement_regions(target_region[["chrom", "start", "end"]])
    # if the target_region do not have all chromosomes, generate full chromosome regions
    done_chrom = set(target_complement["chrom"].tolist())
    empty_chroms = [chrom for chrom in input_df["chrom"].unique().tolist() if chrom not in done_chrom]
    empty_chrom_df = pd.DataFrame({"chrom":empty_chroms,"start":[0]*len(empty_chroms),"end":[1e10]*len(empty_chroms)})
    # finalize the complement regions
    target_complement = pd.concat([target_complement, empty_chrom_df])
    return extract_overlapped_region(input_df, target_complement)

def plot_bed_regions(chrom, start, end, **kwargs):
    temp_df = pd.DataFrame({"chrom":[str(chrom)], "start":[int(start)-10000], "end":[int(end)+10000]})
    df_list = [
        find_overlapped_region(
            input_df=kwargs[k], target_region=temp_df
        )[["chrom", "start", "end"]].assign(region_type=k) for k in kwargs.keys()
    ]
    plot_data = pd.concat(df_list, ignore_index=True)
    chart_list = []
    for i in range(0, plot_data.shape[0]):
        temp_data = plot_data.iloc[[i]].melt(id_vars=["chrom", "region_type"])
        chart = alt.Chart(temp_data).mark_line(size=10).encode(
            x = alt.X("value:Q", scale = alt.Scale(zero=False, domain=[start, end])),
            y = "region_type:N",
            color = "region_type:N"
        )
        chart_list.append(chart)
    chart = chart_list[0]
    for i in range(1, len(chart_list)):
        chart += chart_list[i]
    return chart.interactive()

def get_df_by_region(input_df, chrom, start, end):
    return input_df.loc[
        (input_df["chrom"] == chrom) &
        (input_df["start"] >= start) &
        (input_df["end"] <= end)
    ]