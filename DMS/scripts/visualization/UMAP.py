import os
import pandas as pd
import numpy as np
import pylab as p
from plotnine import *
import matplotlib.pyplot as plt
from scipy.spatial.distance import jensenshannon
from sklearn.metrics import pairwise_distances
import igraph as ig
import umap
import leidenalg as la
import seaborn as sns

#define function
class DMSProfileGraph(ig.Graph):
    def __init__(self, edges_df):
        nodes = pd.unique(edges_df['index'])
        self.sample2node = {}
        self.node2sample = {}
        edges_node = set()
        for i in range(len(nodes)):
            self.node2sample[i] = nodes[i]
            self.sample2node[nodes[i]] = i
        for u, v, w in edges_df.to_numpy():
            _u = self.sample2node[u]
            _v = self.sample2node[v]
            edges_node.add((min(_u, _v), max(_u, _v)))
        super().__init__(n=len(nodes), edges=edges_node)
    def __len__(self):
        return len(self.sample2node)

#load data
dms_results = pd.read_csv("./DMS/files/antibody_dms_merge_avg.csv")
len(dms_results)


site_sum = dms_results.groupby(['site','antibody']).sum().reset_index()
sites_total_sum = site_sum.groupby('antibody')['mut_escape'].sum().to_dict()

site_mat = site_sum.pivot(index='antibody', columns='site', values='mut_escape').fillna(0)
dist = pd.DataFrame(pairwise_distances(site_mat.to_numpy(), metric=jensenshannon), index=site_mat.index, columns=site_mat.index)
dist.index.name = None

n_neighbors = 12
edges = (
    dist.reset_index().melt(id_vars="index")
        .query('index != variable').sort_values('value',ascending=True)
        .groupby('index').head(n_neighbors)
        .reset_index(drop=True)
)
knn_graph = DMSProfileGraph(edges_df = edges)


# UMAP
embedding = pd.DataFrame(
    umap.UMAP(n_neighbors=n_neighbors, min_dist=0.8, metric='precomputed', random_state=42).fit_transform(dist),
    columns=["UMAP1","UMAP2"]
)
embedding['antibody'] = dist.index
embedding.index = [knn_graph.sample2node[i] for i in dist.index]

# Leiden
partition = la.find_partition(knn_graph, partition_type=la.RBConfigurationVertexPartition, resolution_parameter=1, seed=42).membership
embedding['cluster'] = [str(partition[i]) for i in embedding.index]
cluster_counts = embedding['cluster'].value_counts()
print(cluster_counts)

plt.figure(figsize=(4, 4))
sns.scatterplot(data=embedding, x='UMAP1', y='UMAP2', hue='cluster')
plt.show()

cluster_rename = pd.Series(
    [
        'C0',
        'C1',
        'C2',
        'C3'
    ],
    index=np.unique(embedding['cluster'])
).to_dict()

embedding = embedding.sort_values('antibody').reset_index(drop=True).assign(group=lambda x:[cluster_rename[y] for y in x['cluster']])

embedding['group'] = pd.Categorical(embedding['group'], 
                                  categories=['C0', 'C1', 'C2', 'C3'], 
                                  ordered=True)
print(embedding)

highlight_antibodies = ['BRII-196','BD56-1854',
                        'REGN10987','BD57-1303','BD57-1271','BRII-198',
                        'BD55-5514', 'AZD1061']

plt.figure(figsize=(5, 5))
sns.scatterplot(data=embedding, x='UMAP1', y='UMAP2', hue='cluster')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
for i in range(len(embedding)):
    antibody = embedding['antibody'].iloc[i]
    if antibody in highlight_antibodies:
        plt.annotate(
            antibody, 
            xy=(embedding['UMAP1'].iloc[i], embedding['UMAP2'].iloc[i]),
            xytext=(embedding['UMAP1'].iloc[i] + 0.2, embedding['UMAP2'].iloc[i] + 0.2),
            arrowprops=dict(facecolor='black', arrowstyle='->', lw=1),
            fontsize=6, 
            color='black' 
        )
plt.show()

antibody_info = pd.read_csv("./DMS/files/mAb_embed.csv")
mAb_embed = pd.merge(embedding.iloc[:, 0:3], antibody_info, on='antibody', how='left')

mAb_embed.to_csv('./DMS/files/mAb_embed.csv', index=False)
