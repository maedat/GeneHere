#!/usr/bin/env python3
# coding: UTF-8

# 
# ======================================================================
# Project Name    : GeneHere
# File Name       : genehere.py
# Version       : 1.0.0
# Encoding        : python
# Creation Date   : 2019/09/1
# Author : Taro Maeda 
# license     MIT License (http://opensource.org/licenses/mit-license.php)
# Copyright (c) 2019 Taro Maeda
# ======================================================================
# 


#モジュールを読み込む
import argparse
import pandas as pd
from reportlab.lib.units import cm
from Bio.Graphics import BasicChromosome
from Bio import SeqIO
from BCBio import GFF
import random, string
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

def GET_ARGS():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--fasta',  help="File path to a fasta file", required=True)
    parser.add_argument('-g','--gff', help="GFF3 file path", required=True)
    parser.add_argument('-t','--title', help="Figure title (default = No Title)", required=False)
    parser.add_argument('-o','--out', help="Output file (pdf) name (default = simple_chrom.pdf)", required=False)
    
    parser.add_argument('-s','--scale', help="Manualy set the ladder scale size (default = 100000) ", type=int, required=False)
    parser.add_argument('-S','--scale_split', help="How many split ladder (default = 10) ", type=int, required=False)

    parser.add_argument('-x','--xsize', help="output image file size (x axis, default = 100 [cm] )", required=False)
    parser.add_argument('-y','--ysize', help="output image file size (y axis, default = 30 [cm] )", required=False)
    parser.set_defaults(xsize=100, ysize=30)
    parser.set_defaults(title="No Title", out="simple_chrom.pdf")
    parser.set_defaults(scale=100000, scale_split=10)

    return parser.parse_args()
    
    
def RandomSeq(n):
   return ''.join(random.choices("ATGC", k=n))
    

if __name__ == '__main__':

    args = GET_ARGS()
    fasta_in = args.fasta
    gff_in = args.gff
    X_size = args.xsize
    Y_size = args.ysize
    
    title_name = args.title
    output_file = args.out
    
    scale_max = args.scale
    scale_split = args.scale_split
    
    max_len = 0  
    telomere_length = 0  
    
    chr_diagram = BasicChromosome.Organism() 
    chr_diagram.page_size = (X_size*cm, Y_size*cm) 
    df = pd.DataFrame(columns = ['length', 'ID'])
    
    

    with open(fasta_in, "r") as fh:
        for record in SeqIO.parse(fh, "fasta"):
            sdf= pd.DataFrame(
                        {'length': len(record),
                        'ID': record.id},
                        index=[record.id]
                        )
            df=df.append(sdf)
            
    df.length = df.length.astype(int)   
    max_len = df.length.max()
    
    
#scale用ラダー染色体を書く

    lad_rec = SeqRecord(RandomSeq(scale_max), "ladder_"+str(scale_max)+" bp")
    lad_feature = SeqFeature()
    lad_feature_s = []
    
    for i in range(1,scale_split,2):
        region_unit=int(scale_max/scale_split)
        tmp_featur = SeqFeature(FeatureLocation(region_unit*i, region_unit*(i+1)), type="gene", strand=1)
        tmp_start = region_unit*i
        tmp_featur.qualifiers["locus_tag"] = str(tmp_start)
        lad_feature_s.append(tmp_featur)
        
    lad_rec.features = lad_feature_s
    
    cur_chromosome = BasicChromosome.Chromosome(lad_rec.id)
    cur_chromosome.scale_num = max_len + 2 * telomere_length
    
    start = BasicChromosome.TelomereSegment() 
    start.scale = telomere_length
    cur_chromosome.add(start)
    

    feature={}
    body = BasicChromosome.AnnotatedChromosomeSegment(scale_max,lad_rec.features)
    body.scale = scale_max 
    cur_chromosome.add(body)

    end = BasicChromosome.TelomereSegment(inverted=True)
    end.scale = telomere_length
    cur_chromosome.add(end)

    chr_diagram.add(cur_chromosome)


####ladder 終了
####実際の配列スタート


    with open(fasta_in, "r") as fh:
        for record in SeqIO.parse(fh, "fasta"):
            features=[] #featuresを初期化しておく(無い時があるので）
            length = len(record) #配列長を取得する
            ID = record.id #配列名を取得する

            #gffからfeatureを取得する。
            gff_handle = open(gff_in)
            for rec in GFF.parse(gff_handle):
                if rec.id == ID:
                    features = [f for f in rec.features]
                    for f in features:
                        f.qualifiers["color"] = list(map(int, f.qualifiers["color"] ))

            cur_chromosome = BasicChromosome.Chromosome(ID)
            cur_chromosome.scale_num = max_len + 2 * telomere_length 

            # Add an opening telomere
            start = BasicChromosome.TelomereSegment() 
            start.scale = telomere_length
            cur_chromosome.add(start)

            #gffで得たfeature情報を追加する
            # Add a body - using bp as the scale length here.
            body = BasicChromosome.AnnotatedChromosomeSegment(length, features)
            body.scale = length 
            cur_chromosome.add(body)

        # Add a closing telomere
            end = BasicChromosome.TelomereSegment(inverted=True)
            end.scale = telomere_length
            cur_chromosome.add(end)

        # This chromosome is done
            chr_diagram.add(cur_chromosome)


    chr_diagram.draw(output_file, title_name) 
