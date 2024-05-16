# Description: Script (1 out of 2) to analyze the FunC-ESMs classification of the human proteome, the script is divided in steps, each step is a different analysis, the steps are numbered and can be selected with the --step argument
# Script for  the manuscript:"Decoding molecular mechanisms for loss of function variants in the human proteome"
# Author: Matteo Cagiada (matteo.cagiada@bio.ku.dk)
# Date: 2024-05-15
# Version: 1.0
# License: MIT


## run with python 3.9
## list dependencies: numpy, pandas, scipy, vaex, os, sys, argparse, pickle, glob

## import libraries
import numpy as np
import pandas as pd
import scipy
import vaex
import os,sys
import argparse
import pickle
from glob import glob

def create_parser():
    parser = argparse.ArgumentParser(
            prog='human_proteome_analysis_stats.py',
            description='Script to generate the analysis result for the FunC-ESMs classification of the human proteome, the script is divided in steps, each step is a different analysis, the steps are numbered and can be selected with the --step argument'
    )

    parser.add_argument(
        "--step",
        type=int,
        default=0,
        help="Step to run (1-7)"
    )
    parser.add_argument(
        "--out-folder",
        type=str,
        default='./output/',
        help="Output folder"
    )
    parser.add_argument(
        "--df-hdf5",
        type=str,
        default='./functional_proteome.hdf5',
        help="Path to the input dataframe"
    )

    return parser


def main(args):
   
    alphabetAA_L_D={'-':0,'_' :0,'A':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,'K':9,'L':10,'M':11,'N':12,'P':13,'Q':14,'R':15,'S':16,'T':17,'V':18,'W':19,'Y':20}

    ## analysis step
    step=args.step

    # set up the thresholds
    esm1b_threshold=-6.5
    esmif_threshold=-7.0

    # load dataframe from files 
    df_location=args.df_hdf5

    unique_uniprots_path=os.path.join(args.out_folder,'unique_uniprots.npy')

    # create the output folder 
    out_folder=os.path.realpath(args.out_folder)
    
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    

    # load the dataframe
    df_full=vaex.open(df_location)
    
    print('--> Dataframe loaded')
    
    if step == 1:
        # create the uniprot id list from the dataframe
        print('Step 1: Creating unique uniprots')

        unique_uniprots=np.unique(df_full.uniprot.to_numpy())
    
        np.save(os.path.join(out_folder,'unique_uniprots.npy'),unique_uniprots)
        print('--> Unique uniprots saved')

    else:
        unique_uniprots=np.load(unique_uniprots_path,allow_pickle=True)
        print('--> Unique uniprots loaded')
 
    df_res=df_full[df_full['aa_var'] == '=']

    if step ==2:
        # count the variants and residues for the full dataset
        print('Step 2: Counting variants and residues for the full dataset')

        fm_vc=df_full.func_esms_variant_class.to_numpy()
        fm_rc=df_full.func_esms_residue_class.to_numpy()
        aa_type=df_full.aa_var.to_numpy()

        print('-->counting variants')
        counts_vc=np.zeros(4,dtype=int)
        counts_nan=0
        for idx,(aa,vc) in enumerate(zip(aa_type,fm_vc)):
            if aa != '=' and aa!='*':
                if np.isnan(vc):
                    counts_nan+=1
                else:
                    counts_vc[int(vc)]+=1

        print('Variant counts, examinated: ',np.sum(counts_vc), ' out of :', np.sum(counts_vc)+counts_nan )
        print('WT-type: ',counts_vc[0],' - ( ',counts_vc[0]/np.sum(counts_vc),'% )')
        print('SBI: ',counts_vc[1],' - ( ',counts_vc[1]/np.sum(counts_vc),'% )')
        print('total-loss: ',counts_vc[2],' - ( ',counts_vc[2]/np.sum(counts_vc),'% )')
        
        print('-->counting residues')

        counts_rc=np.zeros(6,dtype=int)
        counts_nan=0
        counts_weird=[]
        for idx,(aa,rc) in enumerate(zip(aa_type,fm_rc)):
            if aa == '=':
                if np.isnan(rc):
                    counts_nan+=1
                else:
                    try:
                        counts_rc[int(rc)]+=1
                    except:
                        counts_weird.append(idx)

        print('Residue counts, examinated: ',np.sum(counts_rc), ' out of: ', np.sum(counts_rc)+counts_nan )
        print('Tolerant: ',counts_rc[0],' - ( ',counts_rc[0]/np.sum(counts_rc),'% )')
        print('Functional: ',counts_rc[1],' - ( ',counts_rc[1]/np.sum(counts_rc),'% )')
        print('Structural: ',counts_rc[2],' - ( ',counts_rc[2]/np.sum(counts_rc),'% )')
        print('Mixed: ',counts_rc[4],' - ( ',counts_rc[4]/np.sum(counts_rc),'% )')

    if step ==3:
        # count the variants and residues for the IDR and folded regions
        print('Step 3: Counting variants and residues for the IDR and folded regions')

        count_regions=[0,0,0,0]
        count_regions[0]=len(np.nonzero(df_res['spot_disorder'].to_numpy() == 0.)[0])
        count_regions[1]=len(np.nonzero(df_res['spot_disorder'].to_numpy() == 1.)[0])
        count_regions[2]=len(np.nonzero(np.isnan(df_res['spot_disorder'].to_numpy()))[0])
        count_regions[3]=(len(df_res['spot_disorder'].to_numpy()) - count_regions[0] - count_regions[1]-count_regions[2])
        print('folded regions: ',count_regions[0], 'idr regions: ',count_regions[1], 'nan regions: ', count_regions[2] )
        print(count_regions)
        
        for folded_binary in [0,1]: 
            print('processing class ', folded_binary, 'proteins (1 idr, 0 folded)')

            df_sel=df_full[df_full['spot_disorder'] == folded_binary]
        
            fm_vc=df_sel.func_esms_variant_class.to_numpy()
            fm_rc=df_sel.func_esms_residue_class.to_numpy()
            aa_type=df_sel.aa_var.to_numpy()

            print('-->counting variants')
            counts_vc=np.zeros(4,dtype=int)
            counts_nan=0
            for idx,(aa,vc) in enumerate(zip(aa_type,fm_vc)):
                if aa != '=' and aa!='*':
                    if np.isnan(vc):
                        counts_nan+=1
                    else:
                        counts_vc[int(vc)]+=1

            print('Variant counts, examinated: ',np.sum(counts_vc), ' out of :', np.sum(counts_vc)+counts_nan )
            print('WT-type: ',counts_vc[0],' - ( ',counts_vc[0]/np.sum(counts_vc),'% )')
            print('SBI: ',counts_vc[1],' - ( ',counts_vc[1]/np.sum(counts_vc),'% )')
            print('total-loss: ',counts_vc[2],' - ( ',counts_vc[2]/np.sum(counts_vc),'% )')
            
            print('-->counting residues')

            counts_rc=np.zeros(6,dtype=int)
            counts_nan=0
            counts_weird=[]
            for idx,(aa,rc) in enumerate(zip(aa_type,fm_rc)):
                if aa == '=':
                    if np.isnan(rc):
                        counts_nan+=1
                    else:
                        try:
                            counts_rc[int(rc)]+=1
                        except:
                            counts_weird.append(idx)

            print('Residue counts, examinated: ',np.sum(counts_rc), ' out of: ', np.sum(counts_rc)+counts_nan )
            print('Tolerant: ',counts_rc[0],' - ( ',counts_rc[0]/np.sum(counts_rc),'% )')
            print('Functional: ',counts_rc[1],' - ( ',counts_rc[1]/np.sum(counts_rc),'% )')
            print('Structural: ',counts_rc[2],' - ( ',counts_rc[2]/np.sum(counts_rc),'% )')
            print('Mixed: ',counts_rc[4],' - ( ',counts_rc[4]/np.sum(counts_rc),'% )')

            print('weird:', len(counts_weird),counts_weird[0:10])

    if step == 4:
        ## count amino acid types per class in either idr or folded regions
        print('Step 4: Counting amino acid types per class in either idr or folded regions')

        unique_uniprots=np.load(unique_uniprots_path,allow_pickle=True)
        
        for folded_binary in [0,1]: 
            print('processing class ', folded_binary, 'proteins (1 idr, 0 folded)')
            if folded_binary == 0:
                output_amino_npy=os.path.join(out_folder,'folded_amino_count_classes.npy')
                df_sel=df_full[df_full['spot_disorder'] == folded_binary]
            else:
                output_amino_npy=os.path.join(out_folder,'idr_amino_count_classes.npy')
                df_sel=df_full[df_full['spot_disorder'] == folded_binary]
            
            count_residues=np.zeros((20,5),dtype=int)
            
            amino_acids=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

            for idx,aa in enumerate(amino_acids):
                print('Processing: ',aa)
                df_aa=df_sel[df_sel['aa_ref']==aa]
                residues_class= df_aa[(df_aa['aa_var'] == '=')].func_esms_residue_class.to_numpy()
                print('Found ',len(residues_class),' residues')

                for r in residues_class:
                    if not np.isnan(r):
                        try:
                            count_residues[idx,int(r)]+=1
                        except:
                            print('error: ',idx,r)

            np.save(output_amino_npy,count_residues)

    if step == 5:
        ## count coupled starting and ending amino acid types per class in either idr or folded regions 
        print('Step 5: Counting coupled starting and ending amino acid types per class in either idr or folded regions')
                    
        unique_uniprots=np.load(unique_uniprots_path,allow_pickle=True)

        for folded_binary in [0,1]: 
            print('processing class ', folded_binary, 'proteins (1 idr, 0 folded)')
            if folded_binary == 0:
                output_amino_npy=os.path.join(out_folder,'folded_amino_variant_count_classes.npy')
                df_sel=df_full[df_full['spot_disorder'] == folded_binary]
            else:
                output_amino_npy=os.path.join(out_folder,'idr_amino_variant_count_classes.npy')
                df_sel=df_full[df_full['spot_disorder'] == folded_binary]
            
            count_amino_residues=np.zeros((20,20,6),dtype=int)
            
            amino_acids=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

            for idx,aa in enumerate(amino_acids):
                print('Processing: ',aa)
                df_aa=df_sel[df_sel['aa_ref']==aa]
                
                end_mut=df_aa[(df_aa['aa_var'] != '=') & (df_aa['aa_var'] != '*')].aa_var.to_numpy()
                residues_class= df_aa[(df_aa['aa_var'] != '=') & (df_aa['aa_var'] != '*')].func_esms_variant_class.to_numpy()
                
                for index,(av,r) in enumerate(zip(end_mut,residues_class)):
                    if not np.isnan(r):
                        try:
                            count_amino_residues[idx,int(alphabetAA_L_D[av])-1,int(r)]+=1
                        except:
                            print('error: ',idx,r)

            np.save(output_amino_npy,count_amino_residues)

    if step == 6:
        ## extract exposure values for folded regions from the dataframe
        print('Step 6: Extracting exposure values for folded regions from the dataframe')

        output_csv=os.path.join(out_folder,'folded_stats_exposure_fm_df.csv')
        
        unique_uniprots=np.load(unique_uniprots_path,allow_pickle=True)

        df_folded=df_full[df_full['spot_disorder']==0]
        
        df_folded_res=df_full[(df_folded['spot_disorder']==0) & (df_folded['aa_var'] == '=')]
        
        df_vaex_export=df_folded_res['exposure_SS','exposure_HSE_U','exposure_HSE_D','exposure_rASA','funcModel_predicted_class_residue']

        df_export=df_vaex_export.to_pandas_df()
        df_export.to_csv(output_csv)

    if step ==7:
        ## statistic for cysteines part of disulphide bridges (DSB)
        print('Step 7: Statistic for cysteines part of disulphide bridges (DSB)')


        for cysteine_in_DSB in [0,1] :
            print('processing class ', cysteine_in_DSB, 'cysteine (0 no DBS, 1 yes DSB)')

            df_sel=df_full[(df_full['dsbdetector_is_cys_in_disulfide'] == cysteine_in_DSB) & (df_full['aa_ref'] == 'C')]

            fm_vc=df_sel.func_esms_variant_class.to_numpy()
            fm_rc=df_sel.func_esms_residue_class.to_numpy()
            aa_type=df_sel.aa_var.to_numpy()

            print('-->counting variants')
            counts_vc=np.zeros(4,dtype=int)
            counts_nan=0
            for idx,(aa,vc) in enumerate(zip(aa_type,fm_vc)):
                if aa != '=' and aa!='*':
                    if np.isnan(vc):
                        counts_nan+=1
                    else:
                        counts_vc[int(vc)]+=1

            print('Variant counts, examinated: ',np.sum(counts_vc), ' out of :', np.sum(counts_vc)+counts_nan )
            print('WT-type: ',counts_vc[0],' - ( ',counts_vc[0]/np.sum(counts_vc),'% )')
            print('SBI: ',counts_vc[1],' - ( ',counts_vc[1]/np.sum(counts_vc),'% )')
            print('total-loss: ',counts_vc[2],' - ( ',counts_vc[2]/np.sum(counts_vc),'% )')

            print('-->counting residues')

            counts_rc=np.zeros(6,dtype=int)
            counts_nan=0
            counts_weird=[]
            for idx,(aa,rc) in enumerate(zip(aa_type,fm_rc)):
                if aa == '=':
                    if np.isnan(rc):
                        counts_nan+=1
                    else:
                        try:
                            counts_rc[int(rc)]+=1
                        except:
                            counts_weird.append(idx)

            print('Residue counts, examinated: ',np.sum(counts_rc), ' out of: ', np.sum(counts_rc)+counts_nan )
            print('Tolerant: ',counts_rc[0],' - ( ',counts_rc[0]/np.sum(counts_rc),'% )')
            print('Functional: ',counts_rc[1],' - ( ',counts_rc[1]/np.sum(counts_rc),'% )')
            print('Structural: ',counts_rc[2],' - ( ',counts_rc[2]/np.sum(counts_rc),'% )')
            print('Mixed: ',counts_rc[5],' - ( ',counts_rc[5]/np.sum(counts_rc),'% )')


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    main(args)
