import pandas as pd
import os.path as op
import os
from collections import Counter
import random
import numpy as np

import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import sklearn.metrics

def rarefaction_by_sag(df, color='b', show="genome", tail = 0):
    samples = 0
    clusters = 0

    sampled_clusters = set()
    data = []

    for g, sag in enumerate(random.sample(list(df['sag'].unique()), len(df['sag'].unique()))):
        subdf = df[df['sag'] == sag]

        clusts = subdf['cluster']

        for i, c in enumerate(clusts):
            samples += 1
            if i == (len(clusts) - 1):
                new = 1
            else:
                new = 0
            if c not in sampled_clusters:
                clusters += 1
                sampled_clusters.add(c)

            data.append((samples, clusters, sag, new, g))

    pdf = pd.DataFrame(data, columns=["samples","clusters","sag","sag_start", "sag_number"])
    ones = pdf[pdf['sag_start'] == 1]

    if tail > 10:
        return pdf[pdf['sag_start'] == 1][-tail:]

    if show == "genome":
        #plt.scatter(pdf['samples'],pdf['clusters'], color = color, marker='_')
        #plt.scatter(ones['sag_number'],ones['clusters'], color = color)
        plt.plot('sag_number', 'clusters', data=ones, linestyle='-',marker='o', color = color)

    elif show == "gene":
        plt.scatter(ones['samples'],ones['clusters'], color = color)
        plt.plot(x = ones['samples'],y = ones['clusters'], color = color)

    return pdf



cdir = '/mnt/scgc/simon/simonsproject/gorg-clustering'
identifier = "gorg_sag_orfs"
qual="80minid_m80"
tsv_file = op.join(cdir, "analyses", "{}_{}.tsv".format(identifier, qual))
df = pd.read_csv(tsv_file, sep="\t", names = ['c1','c2'])
df['sag'] = [i.split(".")[0].split("_")[0] for i in df['c2']]

df.rename(columns={'c1':'cluster','c2':'orf'}, inplace=True)

groups = pd.read_excel('/mnt/scgc/simon/simonsproject/info/GORG_16S_basicinfo_20190123.xlsx')
groups.rename(columns={"SAG":'sag'}, inplace=True)

mean_est_group_genome_size = groups.dropna(subset=['Estimated_genome_size']).groupby('Group_short')['Estimated_genome_size'].mean().reset_index().copy()
mean_est_group_genome_size.rename(columns={'Group_short':'Group'}, inplace=True)

df = df.merge(groups[['sag','Group_short']], on='sag', how='left')
df.rename(columns={'Group_short':'Group'}, inplace=True)
df = df.fillna('Unidentified')

from matplotlib.patches import Patch

bigdf2 = pd.DataFrame(columns=["samples","clusters","sag","sag_start", "sag_number", "group","iteration"])

n = 10000
gdict = {}

print("there are {} groups".format(len(df['Group'].unique())))

for i, group in enumerate(df['Group'].unique()):

    subdf = df[df['Group'] == group]
    if group == "Unidentified" or group == "Other" or len(subdf) < n:
        continue
    else:
        for j in range(0, 10):
            print("for {group}, iteration {num}".format(group=group, num=j))
            pdf = rarefaction_by_sag(subdf, show=None, tail=10)
            pdf['group'] = group
            pdf['iteration'] = j
            bigdf2 = pd.concat([bigdf2, pdf])

bigdf2.to_csv("/mnt/scgc/simon/simonsproject/gorg-clustering/analyses/190308_tail_rarefaction_calculations.csv")
