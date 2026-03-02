#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 11:58:12 2023

@author: liuyang
"""

import numpy as np
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




########### load data
ldata_all_concat = sc.read_h5ad('EN_Batch01to08_Velocity_all.h5ad')





########### Brain EndMT selected
adata_Mural_brain_nonGW24_EMT_selected = sc.read_h5ad('EN_brain_Mural_nonGW24_EMT_selected_trajectory.h5ad')

new_adata = sc.AnnData(X=adata_Mural_brain_nonGW24_EMT_selected.raw.X)
new_adata.var_names = adata_Mural_brain_nonGW24_EMT_selected.var_names
adata_Mural_brain_nonGW24_EMT_selected.raw = new_adata

adata_Mural_brain_nonGW24_EMT_selected.layers.update(ldata_all_concat[adata_Mural_brain_nonGW24_EMT_selected.obs_names, adata_Mural_brain_nonGW24_EMT_selected.var_names].layers)
adata_Mural_brain_nonGW24_EMT_selected.obs = pd.concat([adata_Mural_brain_nonGW24_EMT_selected.obs, ldata_all_concat[adata_Mural_brain_nonGW24_EMT_selected.obs_names,].obs.iloc[:, 15:18]], axis=1)




## Preprocess the Data
ldata_Mural_brain = scv.pp.filter_and_normalize(adata_Mural_brain_nonGW24_EMT_selected, min_shared_counts=30, n_top_genes=5000, log=False, copy=True)

scv.pp.neighbors(ldata_Mural_brain, n_neighbors=20, n_pcs=40)
scv.pp.moments(ldata_Mural_brain, n_neighbors=20, n_pcs=40)




## Dynamical Model
scv.tl.recover_dynamics(ldata_Mural_brain, n_jobs=5)
scv.tl.velocity(adata_Mural_brain_nonGW24_EMT_selected_w, mode='dynamical', vkey='velocity_dynamical')
scv.tl.velocity_graph(adata_Mural_brain_nonGW24_EMT_selected_w, vkey='velocity_dynamical', n_jobs=7)
scv.pl.velocity_graph(adata_Mural_brain_nonGW24_EMT_selected_w, vkey='velocity_dynamical', threshold=1, basis='umap_Trajectory', color='newL3_Subtype', save='Brain_nonGW24_EMT_selected/velocity_graph_Brain_Mural_nonGW24_EMT_selected_dynamical.png')




## pseudotime
scv.tl.velocity_pseudotime(adata_Mural_brain_nonGW24_EMT_selected_w, vkey='velocity_dynamical')
scv.pl.scatter(adata_Mural_brain_nonGW24_EMT_selected_w, basis='umap_Trajectory', color='velocity_dynamical_pseudotime', cmap='gnuplot', save='Brain_nonGW24_EMT_selected/velocity_pseudotime_Brain_Mural_nonGW24_EMT_selected_dynamical_original.png')
adata_Mural_brain_nonGW24_EMT_selected_w.obs.drop('root_cells', axis=1, inplace=True)
adata_Mural_brain_nonGW24_EMT_selected_w.obs.drop('end_points', axis=1, inplace=True)



scv.pl.velocity_embedding_stream(adata_Mural_brain_nonGW24_EMT_selected_w, vkey='velocity_dynamical', basis='umap_Trajectory', color='velocity_dynamical_pseudotime', color_map='coolwarm', density=1, title = 'dynamical/stream', save='Brain_nonGW24_EMT_selected/velocity_embedding_Brain_Mural_nonGW24_EMT_selected_dynamical_streamline_density1.png')
scv.pl.velocity_embedding_stream(adata_Mural_brain_nonGW24_EMT_selected_w, vkey='velocity_dynamical', basis='umap_Trajectory', color='velocity_dynamical_pseudotime', color_map='coolwarm', density=2, title = 'dynamical/stream', save='Brain_nonGW24_EMT_selected/velocity_embedding_Brain_Mural_nonGW24_EMT_selected_dynamical_streamline_density2.png')
scv.pl.velocity_embedding_grid(adata_Mural_brain_nonGW24_EMT_selected_w, vkey='velocity_dynamical', basis='umap_Trajectory', color='velocity_dynamical_pseudotime', color_map='coolwarm', arrow_size=1, arrow_length=1, density=1, title = 'dynamical/grid', save='Brain_nonGW24_EMT_selected/velocity_embedding_Brain_Mural_nonGW24_EMT_selected_dynamical_grid_arrowSmall.png')
scv.pl.velocity_embedding_grid(adata_Mural_brain_nonGW24_EMT_selected_w, vkey='velocity_dynamical', basis='umap_Trajectory', color='velocity_dynamical_pseudotime', color_map='coolwarm', arrow_size=1.5, arrow_length=1.5, density=0.8, title = 'dynamical/grid', save='Brain_nonGW24_EMT_selected/velocity_embedding_Brain_Mural_nonGW24_EMT_selected_dynamical_grid_arrowBig.png')




## Identify root cells
adata_Mural_brain_nonGW24_EMT_selected_w.obs['velocity_dynamical_pseudotime_original'] = adata_Mural_brain_nonGW24_EMT_selected_w.obs['velocity_dynamical_pseudotime'].copy()

adata_Mural_brain_nonGW24_EMT_selected_w.obs['Trajectory_CCA_snn_res.1.5'].value_counts()

brain_Mural_nonGW24_EMT_selected_root_cells = adata_Mural_brain_nonGW24_EMT_selected_w.obs[
    (adata_Mural_brain_nonGW24_EMT_selected_w.obs['Trajectory_CCA_snn_res.1.5'] == '6.1') &
    (adata_Mural_brain_nonGW24_EMT_selected_w.obs['newL3_Subtype'] == "pre_EndMT") &
    (adata_Mural_brain_nonGW24_EMT_selected_w.obs['Age'] == "GW06")
].index


for tmp_i, tmp_root_id in enumerate(brain_Mural_nonGW24_EMT_selected_root_cells):
    tmp_root_key = adata_Mural_brain_nonGW24_EMT_selected_w.obs.index.get_loc(tmp_root_id)
    scv.tl.velocity_pseudotime(adata_Mural_brain_nonGW24_EMT_selected_w, vkey='velocity_dynamical', root_key=tmp_root_key)
    tmp_save_path = os.path.join('Mural_root_cells', f'Brain_Mural_nonGW24_EMT_selected_root_{tmp_i:03d}.png')
    scv.pl.scatter(adata_Mural_brain_nonGW24_EMT_selected_w, basis='umap_Trajectory', color='velocity_dynamical_pseudotime', cmap='gnuplot', title=f'root #{tmp_i}  {tmp_root_id}', save=tmp_save_path)



## pseudotime by setting root cells
brain_Mural_nonGW24_EMT_selected_root_cells_key = adata_Mural_brain_nonGW24_EMT_selected_w.obs.index.get_loc(brain_Mural_nonGW24_EMT_selected_root_cells[10])
# change accordingly

scv.tl.velocity_pseudotime(adata_Mural_brain_nonGW24_EMT_selected_w, vkey='velocity_dynamical', root_key=brain_Mural_nonGW24_EMT_selected_root_cells_key)
scv.pl.scatter(adata_Mural_brain_nonGW24_EMT_selected_w, basis='umap_Trajectory', color='velocity_dynamical_pseudotime', cmap='gnuplot', save='Brain_nonGW24_EMT_selected/velocity_pseudotime_Brain_Mural_nonGW24_EMT_selected_dynamical_setRoot.png')


scv.pl.velocity_embedding_stream(adata_Mural_brain_nonGW24_EMT_selected_w, vkey='velocity_dynamical', basis='umap_Trajectory', color='velocity_dynamical_pseudotime', color_map='coolwarm', density=2, title = 'dynamical/stream', save='Brain_nonGW24_EMT_selected/velocity_embedding_Brain_Mural_nonGW24_EMT_selected_dynamical_pseudotime_setRoot_streamline.png')
scv.pl.velocity_embedding_grid(adata_Mural_brain_nonGW24_EMT_selected_w, vkey='velocity_dynamical', basis='umap_Trajectory', color='velocity_dynamical_pseudotime', color_map='coolwarm', arrow_size=1.5, arrow_length=1.5, density=0.8, title = 'dynamical/grid', save='Brain_nonGW24_EMT_selected/velocity_embedding_Brain_Mural_nonGW24_EMT_selected_dynamical_pseudotime_setRoot_grid.png')

scv.pl.velocity_embedding_grid(adata_Mural_brain_nonGW24_EMT_selected_w, vkey='velocity_dynamical', basis='umap_Trajectory', color='velocity_dynamical_pseudotime', color_map='coolwarm', arrow_size=3, arrow_length=2.5, density=0.5, title = '', size=500, figsize=(7, 6), save='Brain_nonGW24_EMT_selected/velocity_embedding_Brain_Mural_nonGW24_EMT_selected_dynamical_pseudotime_setRoot_grid.pdf')





########### Heart EndMT selected
adata_Mural_heart_nonGW24_EMT_selected = sc.read_h5ad('EN_heart_Mural_nonGW24_EMT_selected_trajectory.h5ad')

new_adata = sc.AnnData(X=adata_Mural_heart_nonGW24_EMT_selected.raw.X)
new_adata.var_names = adata_Mural_heart_nonGW24_EMT_selected.var_names
adata_Mural_heart_nonGW24_EMT_selected.raw = new_adata

adata_Mural_heart_nonGW24_EMT_selected.layers.update(ldata_all_concat[adata_Mural_heart_nonGW24_EMT_selected.obs_names, adata_Mural_heart_nonGW24_EMT_selected.var_names].layers)
adata_Mural_heart_nonGW24_EMT_selected.obs = pd.concat([adata_Mural_heart_nonGW24_EMT_selected.obs, ldata_all_concat[adata_Mural_heart_nonGW24_EMT_selected.obs_names,].obs.iloc[:, 15:18]], axis=1)




## Preprocess the Data
ldata_Mural_heart = scv.pp.filter_and_normalize(adata_Mural_heart_nonGW24_EMT_selected, min_shared_counts=30, n_top_genes=5000, log=False, copy=True)

scv.pp.neighbors(ldata_Mural_heart, n_neighbors=20, n_pcs=40)
scv.pp.moments(ldata_Mural_heart, n_neighbors=20, n_pcs=40)



## Dynamical Model
scv.tl.recover_dynamics(ldata_Mural_heart, n_jobs=5)
scv.tl.velocity(adata_Mural_heart_nonGW24_EMT_selected_w, mode='dynamical', vkey='velocity_dynamical')
scv.tl.velocity_graph(adata_Mural_heart_nonGW24_EMT_selected_w, vkey='velocity_dynamical', n_jobs=7)
scv.pl.velocity_graph(adata_Mural_heart_nonGW24_EMT_selected_w, vkey='velocity_dynamical', threshold=1, basis='umap_Trajectory', color='newL3_Subtype', save='Heart_nonGW24_EMT_selected/velocity_graph_Heart_Mural_nonGW24_EMT_selected_dynamical.png')




## pseudotime
scv.tl.velocity_pseudotime(adata_Mural_heart_nonGW24_EMT_selected_w, vkey='velocity_dynamical')
scv.pl.scatter(adata_Mural_heart_nonGW24_EMT_selected_w, basis='umap_Trajectory', color='velocity_dynamical_pseudotime', cmap='gnuplot', save='Heart_nonGW24_EMT_selected/velocity_pseudotime_Heart_Mural_nonGW24_EMT_selected_dynamical_original.png')
adata_Mural_heart_nonGW24_EMT_selected_w.obs.drop('root_cells', axis=1, inplace=True)
adata_Mural_heart_nonGW24_EMT_selected_w.obs.drop('end_points', axis=1, inplace=True)


scv.pl.velocity_embedding_stream(adata_Mural_heart_nonGW24_EMT_selected_w, vkey='velocity_dynamical', basis='umap_Trajectory', color='velocity_dynamical_pseudotime', color_map='coolwarm', density=1, title = 'dynamical/stream', save='Heart_nonGW24_EMT_selected/velocity_embedding_Heart_Mural_nonGW24_EMT_selected_dynamical_streamline_density1.png')
scv.pl.velocity_embedding_stream(adata_Mural_heart_nonGW24_EMT_selected_w, vkey='velocity_dynamical', basis='umap_Trajectory', color='velocity_dynamical_pseudotime', color_map='coolwarm', density=2, title = 'dynamical/stream', save='Heart_nonGW24_EMT_selected/velocity_embedding_Heart_Mural_nonGW24_EMT_selected_dynamical_streamline_density2.png')
scv.pl.velocity_embedding_grid(adata_Mural_heart_nonGW24_EMT_selected_w, vkey='velocity_dynamical', basis='umap_Trajectory', color='velocity_dynamical_pseudotime', color_map='coolwarm', arrow_size=1, arrow_length=1, density=1, title = 'dynamical/grid', save='Heart_nonGW24_EMT_selected/velocity_embedding_Heart_Mural_nonGW24_EMT_selected_dynamical_grid_arrowSmall.png')
scv.pl.velocity_embedding_grid(adata_Mural_heart_nonGW24_EMT_selected_w, vkey='velocity_dynamical', basis='umap_Trajectory', color='velocity_dynamical_pseudotime', color_map='coolwarm', arrow_size=1.5, arrow_length=1.5, density=0.8, title = 'dynamical/grid', save='Heart_nonGW24_EMT_selected/velocity_embedding_Heart_Mural_nonGW24_EMT_selected_dynamical_grid_arrowBig.png')



## Identify root cells
adata_Mural_heart_nonGW24_EMT_selected_w.obs['velocity_dynamical_pseudotime_original'] = adata_Mural_heart_nonGW24_EMT_selected_w.obs['velocity_dynamical_pseudotime'].copy()

adata_Mural_heart_nonGW24_EMT_selected_w.obs['Trajectory_CCA_snn_res.1.8'].value_counts()

heart_Mural_nonGW24_EMT_selected_root_cells = adata_Mural_heart_nonGW24_EMT_selected_w.obs[
    (adata_Mural_heart_nonGW24_EMT_selected_w.obs['Trajectory_CCA_snn_res.1.8'] == '7.0') &
    (adata_Mural_heart_nonGW24_EMT_selected_w.obs['newL3_Subtype'] == "pre_EndMT") &
    (adata_Mural_heart_nonGW24_EMT_selected_w.obs['Age'] == "GW06")
].index


for tmp_i, tmp_root_id in enumerate(heart_Mural_nonGW24_EMT_selected_root_cells):
    tmp_root_key = adata_Mural_heart_nonGW24_EMT_selected_w.obs.index.get_loc(tmp_root_id)
    scv.tl.velocity_pseudotime(adata_Mural_heart_nonGW24_EMT_selected_w, vkey='velocity_dynamical', root_key=tmp_root_key)
    tmp_save_path = os.path.join('Mural_root_cells', f'Heart_Mural_nonGW24_EMT_selected_root_{tmp_i:03d}.png')
    scv.pl.scatter(adata_Mural_heart_nonGW24_EMT_selected_w, basis='umap_Trajectory', color='velocity_dynamical_pseudotime', cmap='gnuplot', title=f'root #{tmp_i}  {tmp_root_id}', save=tmp_save_path)



## pseudotime by setting root cells
heart_Mural_nonGW24_EMT_selected_root_cells_key = adata_Mural_heart_nonGW24_EMT_selected_w.obs.index.get_loc(heart_Mural_nonGW24_EMT_selected_root_cells[2])
# change accordingly

scv.tl.velocity_pseudotime(adata_Mural_heart_nonGW24_EMT_selected_w, vkey='velocity_dynamical', root_key=heart_Mural_nonGW24_EMT_selected_root_cells_key)
scv.pl.scatter(adata_Mural_heart_nonGW24_EMT_selected_w, basis='umap_Trajectory', color='velocity_dynamical_pseudotime', cmap='gnuplot', save='Heart_nonGW24_EMT_selected/velocity_pseudotime_Heart_Mural_nonGW24_EMT_selected_dynamical_setRoot.png')

scv.pl.velocity_embedding_stream(adata_Mural_heart_nonGW24_EMT_selected_w, vkey='velocity_dynamical', basis='umap_Trajectory', color='velocity_dynamical_pseudotime', color_map='coolwarm', density=2, title = 'dynamical/stream', save='Heart_nonGW24_EMT_selected/velocity_embedding_Heart_Mural_nonGW24_EMT_selected_dynamical_pseudotime_setRoot_streamline.png')
scv.pl.velocity_embedding_grid(adata_Mural_heart_nonGW24_EMT_selected_w, vkey='velocity_dynamical', basis='umap_Trajectory', color='velocity_dynamical_pseudotime', color_map='coolwarm', arrow_size=1.5, arrow_length=1.5, density=0.8, title = 'dynamical/grid', save='Heart_nonGW24_EMT_selected/velocity_embedding_Heart_Mural_nonGW24_EMT_selected_dynamical_pseudotime_setRoot_grid.png')

scv.pl.velocity_embedding_grid(adata_Mural_heart_nonGW24_EMT_selected_w, vkey='velocity_dynamical', basis='umap_Trajectory', color='velocity_dynamical_pseudotime', color_map='coolwarm', arrow_size=3, arrow_length=2.5, density=0.5, title = '', size=500, figsize=(7, 6), save='Heart_nonGW24_EMT_selected/velocity_embedding_Heart_Mural_nonGW24_EMT_selected_dynamical_pseudotime_setRoot_grid.pdf')


