#!/usr/bin/env python3
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
print('Importing modules...')
import os
import shutil
import glob
import numpy as np
import pandas as pd
import xlrd as xl
import openpyxl as oxl
from datetime import datetime
print('started ' + datetime.now().strftime("%m/%d/%Y %H:%M:%S"))
#--------------------------------------------------------------------------------------------------#

#------------------------------ VARIABLES ---------------------------------------------------------#
# Not implemented yet...
#--------------------------------------------------------------------------------------------------#

#-------------------------------- FUNCTIONS -----------------------------------------------------------#

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

#------------------------------------------------------------------------------------------------------#

#-------------------------------- EXECUTION -----------------------------------------------------------#
# setting the directory
# if in purdue use below
path = "C:/Users/araujoa/OneDrive/Documentos/Arquivos/Doutorado UESB/work_purdue_1/PhD_Thesis/NSIP_Chapter_4/data_sets/ped_perf"
# if in home use below
#path = "C:/Users/Usuario/Desktop/André/OneDrive/Documentos/Arquivos/Doutorado UESB/work_purdue_1/PhD_Thesis/NSIP_Chapter_4/data_sets/ped_perf"
os.chdir(path) # get in to the folder

#
# Get the data sets to deregress the EBVs
#

print('                  ')
print('Trait 1 running...')
print('                  ')
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


print('                  ')
print('Trait 2 running...')
print('                  ')
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

print('                  ')
print('Trait 3 running...')
print('                  ')
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

print('                  ')
print('Trait 4 running...')
print('                  ')
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

print('                  ')
print('Trait 5 running...')
print('                  ')
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

print('                  ')
print('Trait 6 running...')
print('                  ')
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

print('                  ')
print('Trait 7 running...')
print('                  ')
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


print('                  ')
print('END !!!')
print('Good Luck !!!')
print('finished ' + datetime.now().strftime("%m/%d/%Y %H:%M:%S"))
print('                  ')
####################################################################################################################