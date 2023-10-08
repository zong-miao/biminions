import os
from glob import glob
import pandas as pd
import numpy as np

class probe_regions:
    def __init__(self, probe_bed_file, source = "bed_file"):
        if source == "bed_file":
            self.probe_data = self.fetch_bed_file(probe_bed_file)
            self.probe_data = self.merge_bed_df()
            self.probe_index = self.build_probe_index()
            self.probe_id = [':'.join(data) for data in self.probe_data[["chrom", "start", "end", "GENE"]].astype(str).values.tolist()]
        elif source == "pd_dataframe":
            self.probe_data = probe_bed_file
            self.probe_data = self.merge_bed_df()
            self.probe_index = self.build_probe_index()
            self.probe_id = [':'.join(data) for data in self.probe_data[["chrom", "start", "end", "GENE"]].astype(str).values.tolist()]
        
    def fetch_bed_file(self, file_name):
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
    
    def merge_bed_df(self):
        input_df = self.probe_data
        input_df = input_df.sort_values(["chrom", "start", "end"]).reset_index(drop=True)
        input_df["overlap_group"]=(
            (input_df["start"] >= input_df["end"].shift()) |
            (input_df["chrom"] != input_df["chrom"].shift())
        ).cumsum()
        unique_join = lambda x:';'.join(set(x))
        merged_df = input_df.groupby("overlap_group").agg({
            "chrom":unique_join, "start":"min", "end": "max",
            "GENE":unique_join, "PROBESET": unique_join
        }).reset_index(drop=True)
        return merged_df
    
    def build_probe_index(self):
        """
        self.probe_index has the structure in a dict:
        "chrom:start/10000" : [start, end, i], ...
        """
        posi_index_start = list(self.probe_data["chrom"].astype(str) + ':' +\
        (self.probe_data["start"]/10000).astype(int).astype(str))
        posi_index_end = list(self.probe_data["chrom"].astype(str) + ':' +\
        (self.probe_data["end"]/10000).astype(int).astype(str))
        loci_data = self.probe_data[["start", "end"]].values.tolist()
        
        probe_index = {x:[] for x in list(set(posi_index_start + posi_index_end))}
        for i in range(0, self.probe_data.shape[0]):
            probe_index[posi_index_start[i]].append(loci_data[i] + [i])
            if posi_index_end[i] != posi_index_start[i]:
                probe_index[posi_index_end[i]].append(loci_data[i] + [i])
        return probe_index
    
    def get_probe_index(self, chrom, position):
        result = []
        the_key = f"{chrom}:{int(position/10000)}"
        if the_key not in self.probe_index:
            return result
        for region in self.probe_index[the_key]:
            if position < region[1] and position >= region[0]:
                result.append(region[2])
        return result
    
    def get_probe_id(self, chrom, position):
        chrom = str(chrom)
        position = int(position)
        result = list(set([self.probe_id[i] for i in self.get_probe_index(chrom, position)]))
        return result
    
    def intersect_regions(self, second_regions):
        # merge data #
        combined_df = pd.concat([self.probe_data,second_regions.probe_data]).\
        sort_values(["chrom", "start", "end"]).reset_index(drop = True)
        
        # identify overlapped regions: is this region overlapped with previous row? #
        combined_df["overlapped"] = (
            (combined_df["start"] < combined_df["end"].shift()) &
            (combined_df["chrom"] == combined_df["chrom"].shift())
        )
        
        # get the overlapped regions #
        combined_df["new_start"] = pd.DataFrame({
            "first":combined_df["start"],
            "second":combined_df["start"].shift()
        }).max(axis=1).astype(int)
        combined_df["new_end"] = pd.DataFrame({
            "first":combined_df["end"],
            "second":combined_df["end"].shift()
        }).min(axis=1).astype(int)
        combined_df = combined_df.loc[combined_df["overlapped"] == True]
        
        return combined_df.drop(["start", "end", "overlapped"], axis = 1).rename(
            columns = {"new_start":"start","new_end":"end"}).reset_index(drop=True)   