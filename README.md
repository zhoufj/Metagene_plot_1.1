# Metagene_plot
#Metagene plot for Rec-Seq deep sequencing results. the script was initially provided by Dr. Byung-sik Shin at National Institute of Health for creating metagene plot for visualizing ribosome profiling foot-print distribution across transcripts aligned to main AUG or stop codon. 

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
you need a sam file aligned for transcript, and 
"yeast25_size.txt" having information of utr length
"""

def metagene_len(sam, map_position):
    """get the metagene image plot with 5utr(span5) and 3utr(span3) with
    fragment length information as a heatmap style
    """

    fin1 = open(sam)
    
    Count = 0
    for iR in fin1:
        if iR.startswith('@'):
            Count += 1
        
    dfsam = pd.read_table(sam, sep='\t', header=None, 
                          usecols=[2,3,9],  skiprows=Count)
    dfsam.columns = ['gene','pos','seq']
    dfsam1 = dfsam.assign(Len = dfsam.seq.str.len())
    dfsam2 = dfsam1.loc[:,['gene','pos','Len']]
    dfsam2 = dfsam2.apply(pd.to_numeric, errors='ignore')
    if map_position == 5:
        dfsam4 = dfsam2
    if map_position == 3:
        dfsam3 = dfsam2.assign(pos2 = dfsam2.pos + dfsam2.Len)
        dfsam4 = dfsam3.loc[:,['gene','pos2','Len']]
        dfsam4.columns = ['gene','pos','Len']

    orfsize = pd.read_table('yeast25_size.txt', sep='\t', header=None)
    orfsize.columns = ['gene', 'utr5', 'cds', 'utr3']
      
    dfsam41 = pd.merge(dfsam4, orfsize, on='gene')
    dfsam5 = dfsam41.assign(pos5 = dfsam41.pos - dfsam41.utr5)
    dfsam6 = dfsam5.assign(pos3 = dfsam5.pos - dfsam5.utr5 - dfsam5.cds)
    
    if map_position == 5:
        dfsam7 = dfsam6[(dfsam6['pos5'] >= -15) & (dfsam6['pos5'] <= 5)]
        dfsam8 = dfsam6[(dfsam6['pos3'] <= 10) & (dfsam6['pos3'] >= -10)]
    if map_position == 3:
        dfsam7 = dfsam6[(dfsam6['pos5'] >= -10) & (dfsam6['pos5'] <= 50)]
        dfsam8 = dfsam6[(dfsam6['pos3'] <= 25) & (dfsam6['pos3'] >= -35)]
    dfutr5 = dfsam7.loc[:,['pos5','Len']]
    dfutr3 = dfsam8.loc[:,['pos3','Len']]
    
    def transpose_df(df, Pos):
        df1 = df.groupby([Pos,'Len']).size().reset_index()
        df1.columns = [Pos,'Len','reads']
        df2 = df1[(df1['Len']>=15) & (df1['Len'] <= 34)]
        df3 = pd.pivot_table(df2, values='reads', index=Pos, 
                                 columns='Len', fill_value=0).T
        return df3
    
    dfutr31 = transpose_df(dfutr3, 'pos3')
    dfutr51 = transpose_df(dfutr5, 'pos5')
                             
    return dfutr51, dfutr31

import argparse
import pandas as pd
#import scanpy as sc
import os 



def heatmap_plot(df, name,dirt, vMAX):
                             
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    sns.set(font_scale=2)
    fig, ax = plt.subplots(figsize=(6,3))
    plt.tick_params(axis='x', which='both', top=False)
    plt.tick_params(axis='y', which='both', right=False)
     
    sns.heatmap(df, cmap="Blues", square=True, vmin=0, vmax=vMAX, cbar_kws={"shrink": 0.5}, 
                ax=ax, yticklabels = 5, robust=True,
                linecolor="White").invert_yaxis() 
    plt.savefig(dirt+"/"+name +'.png', dpi=300)

parser = argparse.ArgumentParser(
         description="""
            meta plot.
            """,
            add_help=False
    )

parser.add_argument(
        '-i', '--i_sam',
        action='store',
        dest='i',
        required=True,
        help='input sam file.'
    )

parser.add_argument(
        '-pos', '--pos_num',
        action='store',
        dest='pos',
	    type=int,
        required=True,
        help='postion 5 or 3'
    )

parser.add_argument(
        '-n', '--name_output',
        action='store',
        dest='n',
        help='Layer in anndata to use for plots.'
    )
parser.add_argument(
        '-dirt', '--out_dir',
        action='store',
        dest='dirt',
        help='output directory'
    )

parser.add_argument(
        '-vMAX', '--vMAX_end',
        action='store',
        dest='vMAX',
	    type=int,
        required=True,
        help='The max number of reads in density scale'
    )

options = parser.parse_args()

#if options.output_file != '':
 #       os.makedirs(options.output_file, exist_ok=True)
#
sam=options.i
pos=options.pos
name=options.n
dirt=options.dirt
vmax=options.vMAX
#sam = 'GCATA.sam' ## put your sam file here
df1,df2 = metagene_len(sam, pos) ## 3 for 3' and 5 for 5' end mapping
name5=name + "5_small"
heatmap_plot(df1, name5,dirt,vmax) ## plot for 5' utr. put file name here you want to save
name3=name + "3_small"
heatmap_plot(df2, name3,dirt,vmax) ## plot for 3' utr.
