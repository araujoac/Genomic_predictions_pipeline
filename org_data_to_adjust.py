#!/usr/bin/env python3.9.5
# Name: org_data_to_adjust.py
# Function: Organize the NSIP pre-adujsted data to ajust for comtemporary group (cg)
# Site: See this after
# Author: André C. Araujo
# Manteiner: André C. Araujo
#
#--------------------------------------------------------------------------------------------------#
# Details:  The adjustement of the data is going to be removing the cg mean of preadjusted phenotypes
#           from AGBU. This adjustement is going to be done with R. A function to organize the data
#           for each specific trait was made to do it.
#
#           Before running this program, the input excel sheet sent from Dr Lewis was edited. It was
#           included the class present in the first line in the name of each variable and then the
#           first (class) and second (old header for each trait) were removed. The combination of the 
#           trait name and class was used as the header in the sheete "Performance data" for further
#           analysis.
#
#           At the end, a data set for each trait of interest with the adjobs - cg mean was obtained
#           by using the outputs of this program in the R program called "NSIP_dat_org.R". The columns
#           of the final outputs i this program present: ID, breed code, flock code, cg code, and
#           adjobs phenotypes for each trait. The output from the R program have one colume more beiond
#           the output of this program, that is the adjobs - cg mean.
#
#           Any question fell free to ask Andre C. Araujo: araujoa@pudue.edu
#                                                          andrefuturo7@hotmail.com
#                                                          araujoandrecampelo@gmail.com
#
#--------------------------------------------------------------------------------------------------#
# Change log:
# Criation: v1.0 in 05/12/2021
# Editions: 
#--------------------------------------------------------------------------------------------------#
# Tested in: Python 3.9.5, MSC v.1928 64 bit (AMD64)] on win32
#--------------------------------------------------------------------------------------------------#
# Acknowledges: Brito's Lab - Purdue University
#               GDMA - UVF
#               GACOM - UESB
#--------------------------------------------------------------------------------------------------#

#-------------------------------- IMPORTING MODULES -----------------------------------------------------------#
import os
import shutil
import glob
import numpy as np
import pandas as pd
import xlrd as xl
import openpyxl as oxl
#--------------------------------------------------------------------------------------------------#

#------------------------------ VARIABLES ---------------------------------------------------------#
# Not implemented yet...
#--------------------------------------------------------------------------------------------------#

#-------------------------------- FUNCTIONS -----------------------------------------------------------#

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file
# Parameter b: name of the excel sheet
# Parameter c: name of trait_cg
# Parameter d: name of trait_adjobs
# Parameter e: breed code
# Parameter f: name of the output file     
def prep_perf_data_NSIP(a,b,c,d,e,f):
    print("Reading data...")
    data = pd.read_excel(a, sheet_name=b,
        usecols = ['lpid', 'breed', 'flock', c, d])
    print("Organization ongoing ...")
    breed = data[(data.breed == e)]
    breed_clean = breed.dropna(subset=[c, d])
    breed_ord = breed_clean.sort_values(by=c)
    breed_ord[c] = breed_ord[c].astype(int) # making as int (numeric)
    breed_ord.to_csv(f, index=None, sep=" ", mode="w")
    print(breed_ord.head)
    print("END!!!")
    print("Good Luck!!!")


# Function to prepare the NSIP EBV files to get the dEBVs
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: name of the excel file pedigree data
# Parameter f: name of the excel sheet pedigree data
# Parameter g: breed code
# Parameter h: minimum reliability
# Parameter i: name of the output file     
def prep_data_to_deregress(a,b,c,d,e,f,g,h,i):
    print("Reading EBV data...")
    data_ebv = pd.read_excel(a, sheet_name=b,
        usecols = ['lpid', 'breed', 'flock', c,d,'geno'])
    print("Reading pedigree data...")
    data_ped = pd.read_excel(e, sheet_name=f,
        usecols = ['lpid', 'sire', 'dam', 'breed', 'flock', 'sex'])
    print("Organizing the EBV data...")
    breed_ebv = data_ebv[(data_ebv.breed == g)]
    breed_ebv = breed_ebv.assign(rel=(breed_ebv[d]/100)**2)
    breed_ebv = breed_ebv.dropna(subset=[c,'rel'])
    breed_ebv = breed_ebv[(breed_ebv.rel >= h)]
    print("Organizing the pedigree data...")
    breed_ped = data_ped[(data_ped.breed == g)]
    breed_ped = breed_ped.drop(['breed', 'flock'], axis=1)
    print("Merging the EBV and pedigree data...")
    # using merge function by setting how='inner', commom merge
    breed_ebv_ped = pd.merge(breed_ebv, breed_ped, 
                       on='lpid', 
                       how='inner')             # merging by animal ID (lpid)
    breed_ebv_ped = breed_ebv_ped[['lpid', 'sex', 'sire', 'dam',c, 'rel', 'geno']]
    print("Writing the output...")
    breed_ebv_ped.to_csv(i, index=None, sep=" ", mode="w")
    print("END!!!")
    print("Good Luck!")



# Function to prepare the NSIP dEBV files from CGIL shyne app
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: breed code
# Parameter f: minimum reliability
# Parameter g: name of the output file  
def prep_dEBV_data(a,b,c,d,e,f,g):
    print("Reading the EBV data...")
    data_ebv = pd.read_excel(a, sheet_name=b,
        usecols = ['lpid', 'breed', 'flock', c,d])
    print("Reading the dEBV data...")
    data_debv = pd.read_csv('dEBVs_bwt.csv')
    print("Organizing the EBV data...")
    breed_ebv = data_ebv[(data_ebv.breed == e)]
    breed_ebv = breed_ebv.assign(rel=(breed_ebv[d]/100)**2)
    breed_ebv = breed_ebv.dropna(subset=[c,'rel'])
    breed_ebv = breed_ebv[(breed_ebv.rel >= f)]
    print("Organizing the dEBV data...")
    data_debv = data_debv.rename({'ID': 'lpid'}, axis=1)
    data_ebv_debv = pd.merge(breed_ebv, data_debv, 
                       on='lpid', 
                       how='inner')
    data_ebv_debv = data_ebv_debv[['lpid', 'breed', 'flock', c,'rel','dEBV']]
    data_ebv_debv = data_ebv_debv[~data_ebv_debv.isin([np.nan, np.inf, -np.inf]).any(1)]
    print("Writing the output...")
    data_ebv_debv.to_csv(g, index=None, header=None, sep=" ", mode="w")
    print("END!!!")
    print("Good Luck!")



# Function to prepare the NSIP pedigree data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: breed code
# Parameter d: name of the output file
def prep_ped(a,b,c,d):
    print("Reading pedigree data...")
    data_ped = pd.read_excel(a, sheet_name=b,
        usecols = ['lpid', 'sire', 'dam', 'breed','yob'])
    print("Organizing the pedigree data...")
    breed_ped = data_ped[(data_ped.breed == c)]
    breed_ped = breed_ped.drop(['breed'], axis=1)
    breed_ped = breed_ped.assign(alt_dam = 0)
    breed_ped = breed_ped[['lpid', 'sire', 'dam', 'alt_dam', 'yob']]
    print("Writing the output...")
    breed_ped.to_csv(d, index=None, header=None, sep=" ", mode="w")
    print("END!!!")
    print("Good Luck!")


# Function get the NSIP pedigree data with flock, data birth, and sex information
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: breed code
# Parameter d: name of the output file
def ped_flock_dob(a,b,c,d):
    print("Reading pedigree data...")
    data_ped = pd.read_excel(a, sheet_name=b,
        usecols = ['lpid', 'sire', 'dam', 'breed','flock','dob','sex'])
    print("Organizing the pedigree data...")
    breed_ped = data_ped[(data_ped.breed == c)]
    breed_ped = breed_ped.drop(['breed'], axis=1)
    print("Writing the output...")
    breed_ped.to_csv(d, index=None, header=None, sep=" ", mode="w")
    print("END!!!")
    print("Good Luck!")

# Function to show descritive statistics
# Parameter a: name of the excel file
# Parameter b: name of the excel sheet
# Parameter c: breed code
# Parameter d: name of trait_cg
# Parameter e: name of trait_rawobs
def describe_raw_data(a,b,c,d,e,f):
    print("Reading performance data...")
    raw_data = pd.read_excel(a, sheet_name=b, 
        usecols = ['lpid', 'breed', 'flock', d, e])
    ids_to_use = pd.read_csv(f)
    raw_data = raw_data[(raw_data.breed == c)]
    print(raw_data[e].describe())
    raw_data = raw_data.dropna()
    raw_data = pd.merge(raw_data, ids_to_use, on='lpid', how='inner')
    print("""Descriptive statistics:""")
    print(raw_data[e].describe())





#------------------------------------------------------------------------------------------------------#

#-------------------------------- EXECUTION -----------------------------------------------------------#
# setting the directory
# if in purdue use below
path = "C:/Users/araujoa/OneDrive/Documentos/Arquivos/Doutorado UESB/work_purdue_1/PhD_Thesis/NSIP_Chapter_4/data_sets/ped_perf"
# if in home use below
path = "C:/Users/Usuario/Desktop/André/OneDrive/Documentos/Arquivos/Doutorado UESB/work_purdue_1/PhD_Thesis/NSIP_Chapter_4/data_sets/ped_perf"
os.chdir(path) # get in to the folder

#
# Organize the pedigree files
#

# Function to prepare the NSIP pedigree data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: breed code
# Parameter d: name of the output file
prep_ped('Ramb_ped_perf_data.xlsx','Pedigree data',612,'rambuilet_ped.txt') 

# Function get the NSIP pedigree data with flock, data birth, and sex information
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: breed code
# Parameter d: name of the output file
ped_flock_dob('Ramb_ped_perf_data.xlsx','Pedigree data',612,'ped_flock_dob.txt') 

#
# Check descriptive statistics
#

# Function to show descritive statistics
# Parameter a: name of the excel file
# Parameter b: name of the excel sheet
# Parameter c: breed code
# Parameter d: name of trait_cg
# Parameter e: name of trait_rawobs
describe_raw_data('Ramb_ped_perf_data.xlsx','Performance data',612,'bwt_cg','bwt_rawobs','bwt_ids.txt')
describe_raw_data('Ramb_ped_perf_data.xlsx','Performance data',612,'bwt_cg','bwt_rawobs','bwt_ids_partial.txt')

# Function to show descritive statistics
# Parameter a: name of the excel file
# Parameter b: name of the excel sheet
# Parameter c: breed code
# Parameter d: name of trait_cg
# Parameter e: name of trait_rawobs
describe_raw_data('Ramb_ped_perf_data.xlsx','Performance data',612,'pwt1_cg','pwt1_rawobs','pwt_ids.txt')

# Function to show descritive statistics
# Parameter a: name of the excel file
# Parameter b: name of the excel sheet
# Parameter c: breed code
# Parameter d: name of trait_cg
# Parameter e: name of trait_rawobs
describe_raw_data('Ramb_ped_perf_data.xlsx','Performance data',612,'ywt_cg','ywt_rawobs','ywt_ids.txt')

# Function to show descritive statistics
# Parameter a: name of the excel file
# Parameter b: name of the excel sheet
# Parameter c: breed code
# Parameter d: name of trait_cg
# Parameter e: name of trait_rawobs
describe_raw_data('Ramb_ped_perf_data.xlsx','Performance data',612,'ygfw_cg','ygfw_rawobs','ygfw_ids.txt')

# Function to show descritive statistics
# Parameter a: name of the excel file
# Parameter b: name of the excel sheet
# Parameter c: breed code
# Parameter d: name of trait_cg
# Parameter e: name of trait_rawobs
describe_raw_data('Ramb_ped_perf_data.xlsx','Performance data',612,'yfd_cg','yfd_rawobs','yfd_ids.txt')

# Function to show descritive statistics
# Parameter a: name of the excel file
# Parameter b: name of the excel sheet
# Parameter c: breed code
# Parameter d: name of trait_cg
# Parameter e: name of trait_rawobs
describe_raw_data('Ramb repro perf data.xlsx','Repro data',612,'nlb_cg','nlb','nlb_ids.txt')


#
# Get the data to be adjusted
#

# Traits that will have the data organized: bwt, pwt1, ywt, ygfw, and yfd in the Rambouilet breed (612 code)

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file
# Parameter b: name of the excel sheet
# Parameter c: name of trait_cg
# Parameter d: name of trait_adjobs
# Parameter e: breed code
# Parameter f: name of the output file 
# bwt data here
prep_perf_data_NSIP('Ramb_ped_perf_data.xlsx','Performance data','bwt_cg','bwt_adjobs',612,'bwt_data.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file
# Parameter b: name of the excel sheet
# Parameter c: name of trait_cg
# Parameter d: name of trait_adjobs
# Parameter e: breed code
# Parameter f: name of the output file 
# pwt1 data here
prep_perf_data_NSIP('Ramb_ped_perf_data.xlsx','Performance data','pwt1_cg','pwt1_adjobs',612,'pwt1_data.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file
# Parameter b: name of the excel sheet
# Parameter c: name of trait_cg
# Parameter d: name of trait_adjobs
# Parameter e: breed code
# Parameter f: name of the output file 
# ywt data here
prep_perf_data_NSIP('Ramb_ped_perf_data.xlsx','Performance data','ywt_cg','ywt_adjobs',612,'ywt_data.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file
# Parameter b: name of the excel sheet
# Parameter c: name of trait_cg
# Parameter d: name of trait_adjobs
# Parameter e: breed code
# Parameter f: name of the output file 
# ygfw data here
prep_perf_data_NSIP('Ramb_ped_perf_data.xlsx','Performance data','ygfw_cg','ygfw_adjobs',612,'ygfw_data.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file
# Parameter b: name of the excel sheet
# Parameter c: name of trait_cg
# Parameter d: name of trait_adjobs
# Parameter e: breed code
# Parameter f: name of the output file 
# yfd data here
prep_perf_data_NSIP('Ramb_ped_perf_data.xlsx','Performance data','yfd_cg','yfd_adjobs',612,'yfd_data.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file
# Parameter b: name of the excel sheet
# Parameter c: name of trait_cg
# Parameter d: name of trait_adjobs
# Parameter e: breed code
# Parameter f: name of the output file 
# yfd data here
prep_perf_data_NSIP('Ramb_ped_perf_data.xlsx','Performance data','yfd_cg','yfd_adjobs',612,'yfd_data.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file
# Parameter b: name of the excel sheet
# Parameter c: name of trait_cg
# Parameter d: name of trait_adjobs
# Parameter e: breed code
# Parameter f: name of the output file 
# nlb data here
prep_perf_data_NSIP('Ramb_ped_perf_data.xlsx','Performance data','nlb_cg','nlb_adjobs',612,'nlb_data.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file
# Parameter b: name of the excel sheet
# Parameter c: name of trait_cg
# Parameter d: name of trait_adjobs
# Parameter e: breed code
# Parameter f: name of the output file 
# nlw data here
prep_perf_data_NSIP('Ramb_ped_perf_data.xlsx','Performance data','nlw_cg','nlw_adjobs',612,'nlw_data.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file
# Parameter b: name of the excel sheet
# Parameter c: name of trait_cg
# Parameter d: name of trait_adjobs
# Parameter e: breed code
# Parameter f: name of the output file 
# nlb data here
prep_perf_data_NSIP('Ramb repro perf data.xlsx','Repro data','nlb_cg_new','nlb',612,'nlb_data_pred.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file
# Parameter b: name of the excel sheet
# Parameter c: name of trait_cg
# Parameter d: name of trait_adjobs
# Parameter e: breed code
# Parameter f: name of the output file 
# nlb data here
prep_perf_data_NSIP('Ramb repro perf data.xlsx','Repro data','nlw_cg_new','nlw',612,'nlw_data_pred.txt')

#
# Get the data sets to deregress the EBVs
#

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: name of the excel file pedigree data
# Parameter f: name of the excel sheet pedigree data
# Parameter g: breed code
# Parameter h: minimum reliability
# Parameter i: name of the output file     
prep_data_to_deregress('Ramb EBV accuracies.xlsx','EBV & accuracies','bwt','bwt_acc',
    'Ramb_ped_perf_data.xlsx','Pedigree data',612,0.0625,'bwt_to_deregress.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: name of the excel file pedigree data
# Parameter f: name of the excel sheet pedigree data
# Parameter g: breed code
# Parameter h: minimum reliability
# Parameter i: name of the output file     
prep_data_to_deregress('Ramb EBV accuracies.xlsx','EBV & accuracies','pwt','pwt_acc',
    'Ramb_ped_perf_data.xlsx','Pedigree data',612,0.0625,'pwt_to_deregress.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: name of the excel file pedigree data
# Parameter f: name of the excel sheet pedigree data
# Parameter g: breed code
# Parameter h: minimum reliability
# Parameter i: name of the output file     
prep_data_to_deregress('Ramb EBV accuracies.xlsx','EBV & accuracies','ywt','ywt_acc',
    'Ramb_ped_perf_data.xlsx','Pedigree data',612,0.0625,'ywt_to_deregress.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: name of the excel file pedigree data
# Parameter f: name of the excel sheet pedigree data
# Parameter g: breed code
# Parameter h: minimum reliability
# Parameter i: name of the output file      
prep_data_to_deregress('Ramb EBV accuracies.xlsx','EBV & accuracies','ygfw','ygfw_acc',
    'Ramb_ped_perf_data.xlsx','Pedigree data',612,0.0625,'ygfw_to_deregress.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: name of the excel file pedigree data
# Parameter f: name of the excel sheet pedigree data
# Parameter g: breed code
# Parameter h: minimum reliability
# Parameter i: name of the output file     
prep_data_to_deregress('Ramb EBV accuracies.xlsx','EBV & accuracies','yfd','yfd_acc',
    'Ramb_ped_perf_data.xlsx','Pedigree data',612,0.0625,'yfd_to_deregress.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: name of the excel file pedigree data
# Parameter f: name of the excel sheet pedigree data
# Parameter g: breed code
# Parameter h: minimum reliability
# Parameter i: name of the output file      
prep_data_to_deregress('Ramb EBV accuracies.xlsx','EBV & accuracies','nlb','nlb_acc',
    'Ramb_ped_perf_data.xlsx','Pedigree data',612,0.0625,'nlb_to_deregress.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: name of the excel file pedigree data
# Parameter f: name of the excel sheet pedigree data
# Parameter g: breed code
# Parameter h: minimum reliability
# Parameter i: name of the output file      
prep_data_to_deregress('Ramb EBV accuracies.xlsx','EBV & accuracies','nlw','nlw_acc',
    'Ramb_ped_perf_data.xlsx','Pedigree data',612,0.0625,'nlw_to_deregress.txt')


#
# Merging the EBV with dEBVs data sets (get some important colmuns in the EBV data).
#

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: breed code
# Parameter f: minimum reliability
# Parameter g: name of the output file  
prep_dEBV_data('Ramb EBV accuracies.xlsx','EBV & accuracies','bwt','bwt_acc',612,0.0625,'dEBVs_bwt.txt')


# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: breed code
# Parameter f: minimum reliability
# Parameter g: name of the output file  
prep_dEBV_data('Ramb EBV accuracies.xlsx','EBV & accuracies','pwt','pwt_acc',612,0.0625,'dEBVs_pwt.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: breed code
# Parameter f: minimum reliability
# Parameter g: name of the output file  
prep_dEBV_data('Ramb EBV accuracies.xlsx','EBV & accuracies','ywt','ywt_acc',612,0.0625,'dEBVs_ywt.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: breed code
# Parameter f: minimum reliability
# Parameter g: name of the output file  
prep_dEBV_data('Ramb EBV accuracies.xlsx','EBV & accuracies','ygfw','ygfw_acc',612,0.0625,'dEBVs_ygfw.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: breed code
# Parameter f: minimum reliability
# Parameter g: name of the output file  
prep_dEBV_data('Ramb EBV accuracies.xlsx','EBV & accuracies','yfd','yfd_acc',612,0.0625,'dEBVs_yfd.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: breed code
# Parameter f: minimum reliability
# Parameter g: name of the output file  
prep_dEBV_data('Ramb EBV accuracies.xlsx','EBV & accuracies','nlb','nlb_acc',612,0.0625,'dEBVs_nlb.txt')

# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: breed code
# Parameter f: minimum reliability
# Parameter g: name of the output file  
prep_dEBV_data('Ramb EBV accuracies.xlsx','EBV & accuracies','nlw','nlw_acc',612,0.0625,'dEBVs_nlw.txt')


####################################################################################################################




print("Reading EBV data...")
data_ebv = pd.read_excel('Ramb EBV accuracies.xlsx', sheet_name='EBV & accuracies',
    usecols = ['lpid', 'breed', 'flock', 'bwt','bwt_acc','geno'])
print("Reading pedigree data...")
data_ped = pd.read_excel('Ramb_ped_perf_data.xlsx', sheet_name='Pedigree data',
    usecols = ['lpid', 'sire', 'dam', 'breed', 'flock', 'sex'])

print("Organizing the EBV data...")
breed_ebv = data_ebv[(data_ebv.breed == 612)]
breed_ebv = breed_ebv.assign(rel=(breed_ebv['bwt_acc']/100)**2)
breed_ebv = breed_ebv.dropna(subset=['bwt','rel'])
breed_ebv = breed_ebv[(breed_ebv.rel >= 0)]

print("Organizing the pedigree data...")
breed_ped = data_ped[(data_ped.breed == 612)]
breed_ped = breed_ped.drop(['breed', 'flock'], axis=1)

print("Merging the EBV and pedigree data...")
# using merge function by setting how='inner', commom merge
breed_ebv_ped = pd.merge(breed_ebv, breed_ped, 
                   on='lpid', 
                   how='inner')             # merging by animal ID (lpid)
breed_ebv_ped = breed_ebv_ped[['lpid', 'sex', 'sire', 'dam','bwt', 'rel', 'geno']]

print("Writing the output...")
breed_ebv_ped.to_csv('bwt_to_deregress.txt', index=None, sep=" ", mode="w")





print("Reading the EBV data...")
data_ebv = pd.read_excel('Ramb EBV accuracies.xlsx', sheet_name='EBV & accuracies',
    usecols = ['lpid', 'breed', 'flock', 'bwt','bwt_acc'])
print("Reading the dEBV data...")
data_debv = pd.read_csv('dEBVs_bwt.csv')

print("Organizing the EBV data...")
breed_ebv = data_ebv[(data_ebv.breed == 612)]
breed_ebv = breed_ebv.assign(rel=(breed_ebv['bwt_acc']/100)**2)
breed_ebv = breed_ebv.dropna(subset=['bwt','rel'])
breed_ebv = breed_ebv[(breed_ebv.rel >= 0.0625)]

print("Organizing the dEBV data...")
data_debv = data_debv.rename({'ID': 'lpid'}, axis=1)

data_ebv_debv = pd.merge(breed_ebv, data_debv, 
                   on='lpid', 
                   how='inner')
data_ebv_debv = data_ebv_debv[['lpid', 'breed', 'flock', 'bwt','rel','dEBV']]
data_ebv_debv = data_ebv_debv[(data_ebv_debv.dEBV != 'inf')]
data_ebv_debv = data_ebv_debv[~data_ebv_debv.isin([np.nan, np.inf, -np.inf]).any(1)]

print("Writing the output...")
data_ebv_debv.to_csv('dEBVs_bwt.txt', index=None, header=None, sep=" ", mode="w")
print("END!!!")
print("Good Luck!")



    print("Reading pedigree data...")
    data_ped = pd.read_excel('Ramb_ped_perf_data.xlsx', sheet_name='Pedigree data',
        usecols = ['lpid', 'sire', 'dam', 'breed','yob'])

    print("Organizing the pedigree data...")
    breed_ped = data_ped[(data_ped.breed == 612)]
    breed_ped = breed_ped.drop(['breed'], axis=1)
    breed_ped = breed_ped.assign(alt_dam = 0)
    breed_ped = breed_ped[['lpid', 'sire', 'dam', 'alt_dam', 'yob']]
    breed_ped.to_csv('rambuilet_ped.txt', index=None, header=None, sep=" ", mode="w")