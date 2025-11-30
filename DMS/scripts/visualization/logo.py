import os
import logomaker
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from plotnine import *
from matplotlib.backends.backend_pdf import PdfPages

def plot_res_logo(res, prefix, shownames={}, rownames=None, site_thres=0.0, force_plot_sites=None, force_ylim=None, width=None):
    flat_res = res.pivot(index=['antibody', 'site'], columns='mutation', values='mut_escape').fillna(0)
    sites_total_score = flat_res.sum(axis=1)
    _ = sites_total_score[sites_total_score > site_thres].index
    strong_sites = np.unique(np.array(sorted([i[1] for i in _])))
    print(strong_sites)
    plot_sites = strong_sites
    plot_sites = plot_sites[plot_sites < 520].astype(int)
    print(plot_sites)
    if force_plot_sites is not None:
        plot_sites = force_plot_sites
    flat_res = flat_res[flat_res.index.isin(plot_sites, level=1)]
    _ = pd.DataFrame(sites_total_score)
    _.columns = ['value']
    _['site'] = [i[1] for i in _.index]
    _['antibody'] = [i[0] for i in _.index]
    if rownames is not None:
        Abs = rownames
    else:
        Abs = np.unique([i[0] for i in flat_res.index])
    print(Abs)
    Npages = len(Abs) // 10 + 1
    if width is None:
        width = 30
    with PdfPages(prefix + '_aa_logo.pdf') as pdf:
        for p in range(Npages):
            Abs_p = Abs[p * 10:min(len(Abs), (p + 1) * 10)]
            fig = plt.figure(figsize=(width, len(Abs_p) * 4.6)).subplots_adjust(wspace=0.2, hspace=0.5)
            site2pos = {}
            for i in range(len(plot_sites)):
                site2pos[plot_sites[i]] = i
            for i in range(len(Abs_p)):
                ab = Abs_p[i]
                _ = flat_res.query('antibody == @ab').droplevel(0)
                add_sites = np.setdiff1d(plot_sites, _.index)
                for _site in add_sites:
                    _.loc[_site, :] = 0.0
                _ = _.sort_index()
                _.index = range(len(_))
                ax = plt.subplot(len(Abs_p), 1, i + 1)
                logo = logomaker.Logo(_,
                                      ax=ax,
                                      color_scheme='chemistry',
                                      vpad=.1,
                                      width=.8)
                logo.style_xticks(anchor=1, spacing=1, rotation=90, fontsize=16)
                _max = np.sum(_.to_numpy(), axis=1).max()
                # ax.set_xticklabels(plot_sites[1::2])
                ax.set_xticklabels(plot_sites)
                # ax.set_yticks([])
                ax.tick_params(axis='both', which='both', length=0)
                if force_ylim is not None:
                    ax.set_ylim(0.0, force_ylim)
                ax.yaxis.set_tick_params(labelsize=20)
                if ab in shownames:
                    ax.set_title(shownames[ab], fontsize=16, fontweight="bold")
                else:
                    ax.set_title(ab, fontsize=16, fontweight="bold")
            plt.tight_layout()
            pdf.savefig(bbox_inches='tight',pad_inches=0.5)
            plt.close()

def plot_highlight_res_logo(res, prefix, rownames=None, site_thres=0.0, force_plot_sites=None, force_ylim=None, width=None):
    flat_res = res.pivot(index=['antibody', 'site'], columns='mutation', values='mut_escape').fillna(0)
    sites_total_score = flat_res.sum(axis=1)
    _ = sites_total_score[sites_total_score >= site_thres].index
    strong_sites = np.unique(np.array(sorted([i[1] for i in _])))
    print(strong_sites)
    plot_sites = strong_sites
    plot_sites = plot_sites[plot_sites < 520].astype(int)
    print(plot_sites)
    if force_plot_sites is not None:
        plot_sites = force_plot_sites
    flat_res = flat_res[flat_res.index.isin(plot_sites, level=1)]
    _ = pd.DataFrame(sites_total_score)
    _.columns = ['value']
    _['site'] = [i[1] for i in _.index]
    _['antibody'] = [i[0] for i in _.index]
    if rownames is not None:
        Abs = rownames
    else:
        Abs = np.unique([i[0] for i in flat_res.index])
    print(Abs)
    Npages = len(Abs) // 10 + 1
    if width is None:
        width = 30
    with PdfPages(prefix + '_aa_logo.pdf') as pdf:
        for p in range(Npages):
            Abs_p = Abs[p * 10:min(len(Abs), (p + 1) * 10)]
            fig = plt.figure(figsize=(width, len(Abs_p) * 4.6)).subplots_adjust(wspace=0.2, hspace=0.5)
            site2pos = {}
            for i in range(len(plot_sites)):
                site2pos[plot_sites[i]] = i
            for i in range(len(Abs_p)):
                ab = Abs_p[i]
                _ = flat_res.query('antibody == @ab').droplevel(0)
                add_sites = np.setdiff1d(plot_sites, _.index)
                for _site in add_sites:
                    _.loc[_site, :] = 0.0
                _ = _.sort_index()
                _.index = range(len(_))
                ax = plt.subplot(len(Abs_p), 1, i + 1)
                logo = logomaker.Logo(_,
                                      ax=ax,
                                      color_scheme='chemistry',
                                      vpad=.1,
                                      width=.8)
                logo.style_xticks(anchor=1, spacing=1, rotation=90, fontsize=20) 
                _max = np.sum(_.to_numpy(), axis=1).max()
                
                ax.set_xticklabels(plot_sites, fontsize=20, fontweight='bold')
                
                ax.set_yticklabels([])
                
                ax.tick_params(axis='both', which='both', length=0)
                if force_ylim is not None:
                    ax.set_ylim(0.0, force_ylim)
                
                for spine in ax.spines.values():
                    spine.set_linewidth(2.5)
                
            plt.tight_layout()
            pdf.savefig(bbox_inches='tight',pad_inches=0.5)
            plt.close()

os.chdir("C:/Users/cchan/file/work/cls/Lab/03.WYY_DMS")
df = pd.read_csv("./DMS/files/antibody_dms_merge_bycluster.csv")

df_A = df[df['antibody'].isin(["A"])]
plot_highlight_res_logo(df_A, "./DMS/figures/FigS3D.Mut_Trend/Logo_highlightA", site_thres=2.579376, width=3.4)
df_B = df[df['antibody'].isin(["B"])]
plot_highlight_res_logo(df_B, "./DMS/figures/FigS3D.Mut_Trend/Logo_highlightB", site_thres=2.067766, width=3.4)
df_C = df[df['antibody'].isin(["C"])]
plot_highlight_res_logo(df_C, "./DMS/figures/FigS3D.Mut_Trend/Logo_highlightC", site_thres=3.219874, width=3.4)
df_D1 = df[df['antibody'].isin(["D1"])]
plot_highlight_res_logo(df_D1, "./DMS/figures/FigS3D.Mut_Trend/Logo_highlightD1", site_thres=1.850396, width=3.4)
df_D2 = df[df['antibody'].isin(["D2"])]
plot_highlight_res_logo(df_D2, "./DMS/figures/FigS3D.Mut_Trend/Logo_highlightD2", site_thres=2.499534, width=3.4)
df_E2 = df[df['antibody'].isin(["E2"])]
plot_highlight_res_logo(df_E2, "./DMS/figures/FigS3D.Mut_Trend/Logo_highlightE2", site_thres=2.296488, width=3.4)
df_E3 = df[df['antibody'].isin(["E3"])]
plot_highlight_res_logo(df_E3, "./DMS/figures/FigS3D.Mut_Trend/Logo_highlightE3", site_thres=2.820647, width=3.4)
df_F1 = df[df['antibody'].isin(["F1"])]
plot_highlight_res_logo(df_F1, "./DMS/figures/FigS3D.Mut_Trend/Logo_highlightF1", site_thres=1.948153, width=3.4)
df_F2 = df[df['antibody'].isin(["F2"])]
plot_highlight_res_logo(df_F2, "./DMS/figures/FigS3D.Mut_Trend/Logo_highlightF2", site_thres=3.592498, width=3.4)
df_F3 = df[df['antibody'].isin(["F3"])]
plot_highlight_res_logo(df_F3, "./DMS/figures/FigS3D.Mut_Trend/Logo_highlightF3", site_thres=2.384762, width=3.4)
