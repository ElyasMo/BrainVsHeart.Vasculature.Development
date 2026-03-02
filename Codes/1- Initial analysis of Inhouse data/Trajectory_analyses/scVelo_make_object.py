#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 11:58:12 2023

@author: liuyang
"""

# import numpy as np
import pandas as pd
pd.set_option('display.max_columns', 10)
pd.set_option('display.width', 100)
import anndata as ad
import scanpy as sc
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints# import matplotlib.pyplot as plt
# from matplotlib import rcParams
import scvelo as scv
scv.logging.print_version()
# Running scvelo 0.2.5 (python 3.9.16) on 2023-10-20 10:47.
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

import os





## read .loom file
raw_data_dir = "/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom"
raw_data_files = []
for root, dirs, files in os.walk(raw_data_dir):
    dirs.sort()
    files.sort()
    for file in files:
        if file.endswith(".loom"):
            raw_data_files.append(os.path.join(root, file))

raw_data_files
# ['/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch01/GW24_B.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch01/GW24_H.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch01/GW8_B.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch01/GW9_B.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch02/GW9_H.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch02/Gw6brain.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch02/Gw6heart.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch02/gw11_bra.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch03/gw8b2.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch03/gw8h2.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch04/gw6_Bra_2.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch04/gw7_Bra.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch04/gw7_Hea.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch05/GW12_B.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch05/GW12_H.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch05/GW14_B.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch05/GW14_H.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch05/gw6_Hea_2.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch06/GW10_B2.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch06/GW10_H2.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch06/gw10b.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch06/gw15b.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch08/gw12Bra2.loom',
#  '/home/liuyang/Bioc/EndothelialCellAtlas/Raw_data/scRNA_loom/Batch08/gw12Hea2.loom']


raw_data_loom = {os.path.basename(file_path): scv.read(file_path, cache=True) for file_path in raw_data_files}



raw_data_sampleID = ["B_GW24_03", "H_GW24_03", "B_GW08_01", "B_GW09_02", # Batch01
                     "H_GW09_04", "B_GW06_05", "H_GW06_05", "B_GW11_06", # Batch02
                     "B_GW08_07", "H_GW08_07", # Batch03
                     "B_GW06_09", "B_GW07_08", "H_GW07_08", # Batch04
                     "B_GW12_10", "H_GW12_10", "B_GW14_11", "H_GW14_11", "H_GW06_09", # Batch05
                     "B_GW10_14", "H_GW10_14", "B_GW10_12", "B_GW15_13", # Batch06
                     "B_GW12_15", "H_GW12_15"] # Batch08
raw_data_loom = {new_key: value for new_key, value in zip(raw_data_sampleID, raw_data_loom.values())}




raw_data_nCells = [tmp.n_obs for tmp in raw_data_loom.values()]



adata_all = sc.read_h5ad('EN_Batch01to08_3_analysis.h5ad')




adata_all_cells = adata_all.obs_names
ldata_dict = {}
for i, sampleID in enumerate(raw_data_loom):
    ldata = raw_data_loom[sampleID].copy()
    ldata.obs.index = [sampleID + '_' + index.split(":")[-1] + '-1' for index in ldata.obs.index]
    ldata_dict[sampleID] = ldata[ldata.obs_names.isin(adata_all_cells)]

ldata_nCells = [tmp.n_obs for tmp in ldata_dict.values()]


ldata_dict_filter = {key: value for key, value in ldata_dict.items() if value.n_obs > 0}
ldata_dict_filter = dict(sorted(ldata_dict_filter.items()))


ldata_all_merge_list = []
for key, ldata in ldata_dict_filter.items():
    ldata_all_merge_list.append(scv.utils.merge(adata_all, ldata))


ldata_all_concat = ad.concat(ldata_all_merge_list, axis=0)



ldata_all_concat.write_h5ad('EN_Batch01to08_Velocity_all.h5ad')






