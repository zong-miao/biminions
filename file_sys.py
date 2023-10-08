import os
from glob import glob
import pandas as pd
import pygsheets
import numpy as np
import psycopg2 as pg2
import sys

class fc_file_path:
    def __init__(self, fc_list, fc_dir_list=[
        "/ghds/mrd/data/flowcells/lunar1.3_dev",
        "/ghds/mrd/data/flowcells/lunar1.3_dev/bip313_rc3",
        "/ghds/mrd/data/flowcells/lunar1.3_dev/bip313_rc1_data_freeze",
        "/ghds/mrd/data/flowcells/lunar1.3_dev/bip313_rc1_data_freeze_set2",
        "/ghds/mrd/data/flowcells/lunar1.3_prod/",
        "/ghds/mrd/data/flowcells/lunar1v2_dev",
        "/ghds/ivd/flowcentral"]):
        self.fc_dir_list = fc_dir_list
        # get full path to the flowcell folder #
        self.fc_list = {fcid:{"fc_dir":self.get_fc_dir(fcid)} for fcid in fc_list if len(self.get_fc_dir(fcid)) > 0}
        
        # list all samples in each flowcell #
        for fcid in self.fc_list:
            self.fc_list[fcid]["sample_list"] = self.get_sample_list(fcid)
        
        # list the fcid of each sample #
        self.sample_list = {}
        for fcid in self.fc_list:
            for sample in self.fc_list[fcid]["sample_list"]:
                self.sample_list[sample] = {"fcid":fcid}
                
        # add extra fcid information #
        for fcid in self.fc_list:
            # add month information #
            month_map = {"01": "January", "02": "February", "03": "March", "04": "April",
                        "05": "May", "06": "June", "07": "July", "08": "August",
                        "09": "September", "10": "October", "11": "November", "12": "December"             
            }
            self.fc_list[fcid]["month"] = month_map[fcid[2:4]]
            self.fc_list[fcid]["date"] = '/'.join([fcid[2:4], fcid[4:6]])
            # add probe information #
            probe = "NA"
            for sample_id in self.fc_list[fcid]["sample_list"]:
                temp_file = self.find_file_path(sample_id, ".gh_sample_db.hdr.tsv")
                if os.path.exists(temp_file):
                    probe = (pd.read_csv(temp_file, sep="\t").loc[0, "panel"])
                    break
            self.fc_list[fcid]["panel"] = probe
            
        # add extrac sample information #
        for sample_id in self.sample_list:
            try:
                temp_file = self.find_file_path(sample_id, ".gh_sample_db.hdr.tsv")
                total_n = pd.read_csv(temp_file, sep="\t").loc[0, "on_target_molecule_count"]
                self.sample_list[sample_id]["total_mol_count_M"] = total_n/1e6
            except:
                self.sample_list[sample_id]["total_mol_count_M"] = 0
            
    def get_fc_dir(self, fcid):
        all_file_list = []
        for fc_dir in self.fc_dir_list:
            glob_str = f"{os.path.join(fc_dir, fcid)}*"
            file_list = glob(glob_str)
            file_list = [
                fc_path for fc_path in file_list if (
                    (
                        os.path.isfile(os.path.join(fc_path, "manifest.json")) | 
                        os.path.isfile(os.path.join(fc_path, "autoqc_sample_qc.hdr.tsv"))
                    )
                )
            ]
            all_file_list.extend(file_list)
        if len(all_file_list) > 0:
            all_file_list.sort(
                key=lambda fc_path: os.path.getmtime(
                    os.path.join(fc_path, "autoqc_sample_qc.hdr.tsv")
                )
            )
            return all_file_list[-1]
        return ''
    
    def get_sample_list(self, fcid):
        if fcid not in self.fc_list:
            return []
        fc_dir = self.fc_list[fcid]["fc_dir"]
        samples = [x.split('/')[-1].split('.')[0] for x in glob(f"{fc_dir}/*.annotated.somatic.tsv")]
        return list(set(samples))

    def find_file_path(self, sample_id, file_suffix):
        if sample_id not in self.sample_list:
            return ""
        fc_dir = self.fc_list[self.sample_list[sample_id]["fcid"]]["fc_dir"]
        glob_str = f"{fc_dir}/{sample_id}{file_suffix}"
        file_path_list = glob(glob_str)
        if len(file_path_list):
            return file_path_list[0]
        else:
            return ""
    
    def get_fcid(self, sample_id):
        return self.sample_list[sample_id]["fcid"]
    
    def get_fc_trait(self, fcid, target_trait):
        return self.fc_list[fcid][target_trait]
    
    def get_sample_trait(self, sample_id, target_trait):
        return self.sample_list[sample_id][target_trait]

# Fetch db credentials
sys.path.append("/ghds/groups/lunar/notebooks/jkurata/gh_sample_tracker/sampletracker/sampletracker/")
from env import SAMP_DB_USERNAME, SAMP_DB_PASSWORD, TEST_DB_USERNAME, TEST_DB_PASSWORD
def lunardb_sql_query_to_df(query):
    """Run an SQL query on the Guardant PostgreSQL server and return the
    result as a pandas dataframe.
    Args:
        query (str): SQL query string.
    Returns:
        pandas.core.frame.DataFrame: pandas dataframe containing SQL
        query.
    """
    connect = pg2.connect(host="ghbi-live-lunar-sample-metadata-db.clrdmmintk6b.us-west-2.rds.amazonaws.com",
                          database="lunar_sample_metadata_db",
                          user=SAMP_DB_USERNAME,
                          password=SAMP_DB_PASSWORD)
    df = pd.read_sql_query(query, connect)
    connect.close()
    return df

def ghdb_sql_query_to_df(query):
    """Run an SQL query on the Guardant PostgreSQL server and return the
    result as a pandas dataframe.
    Args:
        dbinfo (namedtuple): Namedtuple containing postgres database
                             info (see below).
        query (str): SQL query string.
    Returns:
        pandas.core.frame.DataFrame: pandas dataframe containing SQL
        query.
    """
    connect = pg2.connect(host="10.4.170.26",
                     database="ghdb",
                     user="admin",
                     password="N7Tks0xPS")
    df = pd.read_sql_query(query, connect)
    connect.close()
    return df

def gsheet_to_df(gsheet, worksheet, start_pos="A1"):
    """Use the GoogleDocs python API to obtain a specific 
    worksheet within a GoogleSheet and return it as a pandas
    DataFrame.

    Args:
        gsheet (str): Name of Google sheet in Guardant gsheets.
        worksheet (str): Work sheet (i.e. tab to download as df)

    Returns:
        pandas.core.frame.DataFrame: pandas dataframe containing the worksheet.
    """
    gc = pygsheets.authorize(credentials_directory = "/home/zmiao/MyKeys")
    sh = gc.open(gsheet)
    wsh = sh.worksheet_by_title(worksheet)
    df = wsh.get_as_df(empty_value=np.nan, start=start_pos, numerize=False)
    return df

def df_upload_to_gsheet(input_df, gsheet, worksheet, start_pos="A1"):
    """Use the GoogleDocs python API to upload a data frame
    to a specific worksheet within a GoogleSheet
    Args:
        gsheet (str): Name of Google sheet in Guardant gsheets.
        worksheet (str): Work sheet (i.e. tab to upload to)

    Returns:
        None
    """
    gc = pygsheets.authorize(credentials_directory = "/home/zmiao/MyKeys")
    sh = gc.open(gsheet)
    wsh = sh.add_worksheet(worksheet)
    wsh.set_dataframe(input_df, start=start_pos)
    return