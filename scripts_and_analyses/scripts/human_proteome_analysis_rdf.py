# Description: Script (2 out of 2) to analyze the FunC-ESMs classification of the human proteome, the script is divided in steps, each step is a different analysis, the steps are numbered and can be selected with the --step argument
# Script for the manuscript:"Decoding molecular mechanisms for loss of function variants in the human proteome"
# Author: Matteo Cagiada (matteo.cagiada@bio.ku.dk)
# Date: 2024-05-15
# Version: 1.0
# License: MIT

##run with python 3,9
## list dependencies: numpy, pandas, scipy, vaex, os, sys, argparse, pickle, mdtraj, MDAnalysis, glob

## import libraries
import numpy as np
import pandas as pd
import scipy
import vaex
import os,sys
import argparse
import pickle
import mdtraj as md
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from glob import glob
import warnings

def create_parser():
    parser = argparse.ArgumentParser( prog="human_proteome_analysis_rdf.py",
                                     description='Script to generate clustering analysis, including radial distribution function analysis, for the FunC-ESMs classification of the human proteome, the script is divided in steps, each step is a different analysis, the steps are numbered and can be selected with the --step argument'
    )

    parser.add_argument(
        "--step",
        type=int,
        default=0,
        help="Step to run (1-3"
    )
    parser.add_argument(
        "--out-folder",
        type=str,
        default='./output/',
        help="Output folder"
    )
    parser.add_argument(
        "--df-hdf5",
        default='./functional_proteome.hdf5',
        help="Path to the input dataframe"
    )
    parser.add_argument(
        "--af2-folder",
        default='./alphafold_v2/',
        help="Path to the alphafold2 PDB folder"
    )

    return parser


## standard deviation and mean cummulative classes

class SD_cum(object):
    ## class to evaluate standard deviation in a cumulative way
    def __init__(self,bins,n_classes):
        self.sum=np.zeros((bins,n_classes),dtype=float)
        self.sum2=np.zeros((bins,n_classes),dtype=float)
        self.n=np.zeros((bins,n_classes),dtype=float)
        self.bins=bins
        self.n_classes=n_classes
        
    def eval(self, x, c):
        sd_array=[]
        for b in range(self.bins):
            if self.n_classes == 1:
                self.sum[b]  += x[b]
                self.sum2[b] += x[b]*x[b]
                self.n[b]    += 1.0
                sum, sum2, n = self.sum[b], self.sum2[b], self.n[b]            
            else:
                self.sum[b,c]  += x[b]
                self.sum2[b,c] += x[b]*x[b]
                self.n[b,c]    += 1.0
                sum, sum2, n = self.sum[b,c], self.sum2[b,c], self.n[b,c]
            sd_array.append(np.sqrt(sum2/n - sum*sum/n/n))
        return np.array(sd_array)

class MEAN_cum(object):
    ## class to evaluate mean in a cumulative way
    def __init__(self,bins,n_classes):
        self.sum=np.zeros((bins,n_classes),dtype=float)
        self.n=np.zeros((bins,n_classes),dtype=float)
        self.bins=bins
        self.n_classes=n_classes
    def eval(self, x,c):
        mean_array=[]
        for b in range(self.bins):
            if self.n_classes == 1:
                self.sum[b]  += x[b]
                self.n[b]    += 1.0
                sum, n = self.sum[b], self.n[b]
            else:
                self.sum[b,c]  += x[b]
                self.n[b,c]    += 1.0
                sum, n = self.sum[b,c], self.n[b,c]
            mean_array.append(sum/n)
        return np.array(mean_array)
    
    def return_n(self):
        return self.n


def protein_size_class(in_array):
    ## classifies protein based on protein length

    ## input arra: size 3 array it includes: len_seq, max dist, number_funct sites)
    ## 0-100 small proteins
    ## 100 -250 single-domain-like proteins
    ## 250-500 multi-domains proteins
    ### > 500 xlarge multi-domains protein:
    
    len_prot=in_array[0]
    if len_prot < 100:
        return 0
    elif len_prot <250:
        return 1
    elif len_prot <500:
        return 2
    else:
        return 3

def mask_folded(pos_list,len_protein):
    ## mask the distance matrix based on the pos list given as input

    out_mask=np.ones((len_protein,len_protein),dtype=bool)
    
    ##if missmatch between uniprot lenght and af2 just skip the protein masking everything.
    try:
        for idx1,aa1 in enumerate(pos_list):
            for idx2,aa2 in enumerate(pos_list[(idx1+1):]):
                out_mask[aa1-1,aa2-1]=0
    except:
        out_mask=np.ones((len_protein,len_protein),dtype=bool)
    
    return out_mask

def evaluate_distances(pdb_loc,mask_folded=None):
    ## evaluate COM distance matrix from a pdb file 
    u = mda.Universe(pdb_loc, pdb_loc)

    # Calculate COM positions of residues
    com_positions = []
    for residue in u.residues:
        com_position = residue.atoms.center_of_mass()
        com_positions.append(com_position)
    com_positions=np.array(com_positions)
    # Calculate distance matrix
    distances = distance_array(com_positions, com_positions)
    return distances

def Rg(pdb_loc,res_idx):
    ## evaluate radius of gyration of a protein based on a pdb file and a list of residue index
    try:
        pdb=md.load(pdb_loc)
    
        list_residx = output = [str(x-1) for x in res_idx]
        str_idx=" ".join(list_residx)
    
        topology=pdb.topology
        chainCA=topology.select('chainid 0 and protein and name CA and resid '+str_idx)

        pdb_CA=pdb.atom_slice(chainCA)
   
        rg=md.compute_rg(pdb_CA)[0]*10
    except:
        rg =-99999
    return rg

def main(args):
    warnings.filterwarnings("ignore")    
    
    ## set the step to run
    step=args.step

    ## load the dataframe
    df_location=os.path.abspath(args.df_hdf5)

    af2_folder=os.path.abspath(args.af2_folder)
    
    variants_class_column='/projects/prism/people/bqm193/runs/proteome_analysis/esms_simple_v1/variant_class_threshold_model.npy'
    residues_class_column='/projects/prism/people/bqm193/runs/proteome_analysis/esms_simple_v1/res_class_threshold_model.npy'
    #histogram params
    
    out_folder=os.path.realpath(args.out_folder)
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    ## load the dataframe
    df_full=vaex.open(df_location)

    print('--> Dataframe loaded')

    print('removing idrs...')

    df_aa_folded=df_full[(df_full['spot_disorder']==0) & (df_full['aa_var']=='=')]

    df_full_wt=df_full[(df_full['aa_var']=='=')]

    if step==0:
        # create the uniprot id list of folded proteins (proteins with over 50% folded residues by SPOT)

        unique_uniprots=np.unique(df_full_wt.uniprot.to_numpy())
        folded_proteins=[]

        for idx,uniprot in enumerate(unique_uniprots):
            
            if (idx+1)%100 ==0:
                print(idx+1,' out of ',len(unique_uniprots), ' processed' )

            len_full=len(df_full[(df_full_wt['uniprot']==uniprot)])
            len_folded=len(df_full[(df_full_wt['uniprot']==uniprot) & (df_full_wt['spot_disorder']==0)])

            if len_folded/len_full> 0.5:
                folded_proteins.append(uniprot)

        np.save(os.path.join(out_folder,'unique_uniprots_folded.npy'),folded_proteins)
        print('--> Unique uniprots saved')
    
    if step ==1:
        ## evaluate pairwise distances for each protein in the list of folded proteins

        unique_uniprots=np.load(os.path.join(out_folder,'unique_uniprots_folded.npy'),allow_pickle=True)

        out_npz_folder=os.path.join(out_folder,'pairwise_output_rg_small_bins')
        if not os.path.exists(out_npz_folder):
            os.mkdir(out_npz_folder)

        ## loop over the all the folded proteins, you can split the list in 4 parts and run multiple scripts in parallel to speed up the process
        for idx,uniprot in enumerate(unique_uniprots[:]):
        
            if (idx+1)%10 ==0:
                print(idx+1,' out of ',len(unique_uniprots), ' processed' )
            
            output_collected_npz=os.path.join(out_npz_folder,(uniprot+'_pairwise_distrs.npz'))

            if os.path.exists(output_collected_npz):
                print(uniprot,': evaluation already processed')
                continue

            df_uniprot_aa_folded=df_aa_folded[(df_aa_folded['uniprot']==uniprot)].to_pandas_df()
        
            if glob(os.path.join(os.path.abspath(af2_folder),'*-'+uniprot+'-*'+'.pdb')) ==[]:
                print(uniprot, ' - no pdb found')
                continue
            else:
                af2_struct=glob(os.path.join(os.path.abspath(af2_folder),'*'+uniprot+'*'+'.pdb'))[0]
            
            distmap=evaluate_distances(af2_struct)
           
            ## functionally relevant distances
            
            df_sel=df_uniprot_aa_folded[df_uniprot_aa_folded['func_esms_residue_class']==1.]
            
            number_funct=len(df_sel)
            
            pos_list=df_sel.resi.tolist()

            print(uniprot,' functionally relevant positions: ',[i+1 for i in pos_list])
            
            sel_mask=mask_folded(pos_list,distmap.shape[0])
            
            masked_array=np.ma.masked_array(distmap, mask=sel_mask)

            real_funct_pairwise= masked_array[~masked_array.mask]
            
            ## structurally critical distances
            
            df_sel=df_uniprot_aa_folded[df_uniprot_aa_folded['func_esms_residue_class']==2.]
            
            number_struct=len(df_sel)
            
            pos_list=df_sel.resi.tolist()
            print(uniprot,' tl pos: ',[i+1 for i in pos_list])

            sel_mask=mask_folded(pos_list,distmap.shape[0])

            masked_array=np.ma.masked_array(distmap, mask=sel_mask)
            
            real_struct_pairwise= masked_array[~masked_array.mask]
            
            ## shuffled functional site distances
           
            df_uniprot_aa_folded_notl=df_uniprot_aa_folded[df_uniprot_aa_folded['func_esms_variant_class'] != 2.]
            df_uniprot_aa_folded_notl['funcModel_predicted_class_residue_shuffled']=np.random.permutation(df_uniprot_aa_folded_notl['func_esms_residue_class'].values)

            df_sel=df_uniprot_aa_folded_notl[df_uniprot_aa_folded_notl['funcModel_predicted_class_residue_shuffled']==1.]
            
            pos_list=df_sel.resi.tolist()

            sel_mask=mask_folded(pos_list,distmap.shape[0])

            masked_array=np.ma.masked_array(distmap, mask=sel_mask)

            shuffled_funct_pairwise= masked_array[~masked_array.mask]

            ## evalute rg
            all_pos=df_uniprot_aa_folded.resi.tolist()
            rg_prot=Rg(af2_struct,all_pos)
            
            ## create extra info array
            ### in order:  length protein,  max pairwise distance,, number of functional sites,number of tl sites
            extra_info=[distmap.shape[0],np.nanmax(distmap),number_funct,number_struct,rg_prot]
            
            ## save npz with all the data
            np.savez(output_collected_npz,
                    real_funct_pairwise=real_funct_pairwise,
                    real_struct_pairwise=real_struct_pairwise,
                    shuffled_funct_pairwise=shuffled_funct_pairwise,
                    info_protein=np.array(extra_info))

    if step==2:
        ## evaluate the protein size and the number of functional and total-loss sites for each protein in the list of folded proteins

        in_npz_folder=os.path.join(out_folder,'pairwise_output_rg_small_bins')
        unique_uniprots=np.load(os.path.join(out_folder,'unique_uniprots_folded.npy'),allow_pickle=True)

        uniprot_analysed=[]
        protein_size=[]
        max_internal_dist=[]
        number_funct_sites=[]
        number_struct_sites=[]

        for idx,uniprot in enumerate(unique_uniprots[:]):
        
            if (idx+1)%100 ==0:
                print(idx+1,' out of ',len(unique_uniprots), ' processed' )

            if glob(os.path.join(os.path.abspath(in_npz_folder),uniprot+'*'+'.npz')) ==[]:
                print(uniprot, ' - no npz found')
                continue     
            else:
                in_data_loc=glob(os.path.join(os.path.abspath(in_npz_folder),'*'+uniprot+'*'+'.npz'))[0] 

            ## reading the npz file generated in step 1
            pairwise_data=np.load(in_data_loc)

            uniprot_analysed.append(uniprot)
            protein_size.append(pairwise_data['info_protein'][0])
            max_internal_dist.append(pairwise_data['info_protein'][1])
            number_funct_sites.append(pairwise_data['info_protein'][2])
            number_struct_sites.append(pairwise_data['info_protein'][3])
            
            out_npz=os.path.join(out_folder,'statistic_funct_protein_size.npz')
            np.savez(out_npz,
                    uniprot_analyzed=np.array(uniprot_analysed),
                    protein_size=np.array(protein_size),
                    max_internal_dist=np.array(max_internal_dist),
                    number_funct_sites=np.array(number_funct_sites),
                    number_struct_sites=np.array(number_struct_sites))

    if step ==3:
        ## evaluate the radial distribution mean  function for the functional and total-loss sites for each protein in the list of folded proteins
        
        in_npz_folder=os.path.join(out_folder,'pairwise_output_rg_small_bins')
        unique_uniprots=np.load(os.path.join(out_folder,'unique_uniprots_folded.npy'),allow_pickle=True)

        #bins_features
        bins_range=np.arange(0,30.1,0.1)
        bins_mean=(bins_range[:-1]+bins_range[1:])/2
        bins_thickness=(bins_range[1:]-bins_range[:-1])
        
        #normalization factor
        volume_bin=4*np.pi*np.multiply(np.square(bins_mean),bins_thickness)
        
        # 4 classes of protein based on seq lenght for now (maybe on max distance later)

        real_funct_sites_hist_cum=np.zeros((len(bins_mean),4),dtype=float)
        real_funct_sites_hist_mean=np.zeros((len(bins_mean),4),dtype=float)
        real_funct_sites_hist_std=np.zeros((len(bins_mean),4),dtype=float)
        shuffled_funct_sites_hist_cum=np.zeros((len(bins_mean),4),dtype=float)
        shuffled_funct_sites_hist_mean=np.zeros((len(bins_mean),4),dtype=float)
        shuffled_funct_sites_hist_std=np.zeros((len(bins_mean),4),dtype=float)
        real_struct_sites_hist_cum=np.zeros((len(bins_mean),4),dtype=float)
        real_struct_sites_hist_mean=np.zeros((len(bins_mean),4),dtype=float)
        real_struct_sites_hist_std=np.zeros((len(bins_mean),4),dtype=float)

        mean_fr_cum = MEAN_cum(len(bins_mean),4)
        mean_sc_cum = MEAN_cum(len(bins_mean),4)
        mean_cg_cum = MEAN_cum(len(bins_mean),4)

        sd_fr_cum = SD_cum(len(bins_mean),4)
        sd_sc_cum = SD_cum(len(bins_mean),4)
        sd_cg_cum = SD_cum(len(bins_mean),4)

        count_class=[0,0,0,0]

        for idx,uniprot in enumerate(unique_uniprots[:]):
            if (idx+1)%10 ==0:
                print(idx+1,' out of ',len(unique_uniprots), ' processed' )
            
            if glob(os.path.join(os.path.abspath(in_npz_folder),uniprot+'*'+'.npz')) ==[]:
                print(uniprot, ' - no npz found')
                continue
            else:
                in_data_loc=glob(os.path.join(os.path.abspath(in_npz_folder),'*'+uniprot+'*'+'.npz'))[0]

            pairwise_data=np.load(in_data_loc)
            
            if pairwise_data['info_protein'][2] <=10:
                print('not enough functional residues')
                continue 
            if pairwise_data['info_protein'][3] ==0:
                print('no total-loss residues')
                continue 
            if pairwise_data['info_protein'][4] <0:
                print('rg negative')
                continue 
            
            class_prot=protein_size_class(pairwise_data['info_protein'])
            
            count_class[class_prot]+=1
             
            density_prot=float(pairwise_data['info_protein'][2])/(4/3*np.pi*np.power(pairwise_data['info_protein'][4],3))

            ## checks for possible nan values in the pairwise data

            if np.isnan(pairwise_data['real_funct_pairwise']).any():
                print(uniprot,' - nan in real funct')
                continue

            #real funct hist
            
            counts,bins=np.histogram(pairwise_data['real_funct_pairwise'],bins=bins_range)
            
            ## skip the protein if found without functional residues
            if np.sum(counts) ==0:
                count_class[class_prot]-=1
                print(uniprot,' - no functional in range')
                continue
            
            counts_N=counts/np.square(float(pairwise_data['info_protein'][2]))*2 # time 2 because I am evaluating on half of the distance matrix, but I should count for each position all tneighbour, basically using a complete row
            
            #normalization by volume and density
            counts_Nvol=counts_N/(volume_bin*density_prot)
            
            real_funct_sites_hist_mean[:,class_prot]=mean_fr_cum.eval(counts_Nvol,class_prot)
            real_funct_sites_hist_std[:,class_prot]=sd_fr_cum.eval(counts_Nvol,class_prot)

            #shuffle funct hist
            
            counts,bins=np.histogram(pairwise_data['shuffled_funct_pairwise'],bins=bins_range)
            
            counts_N=counts/np.square(float(pairwise_data['info_protein'][2]))*2 # time 2 because I am evaluating on half of the distance matrix, but I should count for each position all tneighbour, basically using a complete row
           
            #normalization by volume and density
            counts_Nvol=counts_N/(volume_bin*density_prot)

            shuffled_funct_sites_hist_mean[:,class_prot]=mean_cg_cum.eval(counts_Nvol,class_prot)
            shuffled_funct_sites_hist_std[:,class_prot]=sd_cg_cum.eval(counts_Nvol,class_prot)

            #real structural hist

            counts,bins=np.histogram(pairwise_data['real_struct_pairwise'],bins=bins_range)
            
            counts_N=counts/np.square(float(pairwise_data['info_protein'][3]))*2 # time 2 because I am evaluating on half of the distance matrix, but I should count for each position all tneighbour, basically using a complete row

            #normalization by volume and density
            counts_Nvol=counts_N/(volume_bin*density_prot)
            
            real_struct_sites_hist_mean[:,class_prot]=mean_sc_cum.eval(counts_Nvol,class_prot)
            real_struct_sites_hist_std[:,class_prot]=sd_sc_cum.eval(counts_Nvol,class_prot)

        out_hist=os.path.join(out_folder,'hist_cluster_analysis_human_proteome_rg.npz')
        
        np.savez(out_hist,
                real_funct_sites_hist_mean=real_funct_sites_hist_mean,
                real_funct_sites_hist_std=real_funct_sites_hist_std,
                shuffled_funct_sites_hist_mean=shuffled_funct_sites_hist_mean,
                shuffled_funct_sites_hist_std=shuffled_funct_sites_hist_std,
                real_struct_sites_hist_mean=real_struct_sites_hist_mean,
                real_struct_sites_hist_std=real_struct_sites_hist_std,
                bins=bins_range)

        ### save data
        print('total protein manual: ',count_class) 

if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    main(args)
