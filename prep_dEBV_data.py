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
# Parameter c: name of the excel sheet ebv data
# Parameter d: name of trait_ebv
# Parameter e: name of trait_acc
# Parameter f: breed code
# Parameter g: minimum reliability
# Parameter h: name of the output file  
def prep_dEBV_data(a,b,c,d,e,f,g,h):
    print("Reading the EBV data...")
    data_ebv = pd.read_excel(a, sheet_name=b,
        usecols = ['lpid', 'breed', 'flock', d,e])
    print("Reading the dEBV data...")
    data_debv = pd.read_csv(c)
    print("Organizing the EBV data...")
    breed_ebv = data_ebv[(data_ebv.breed == f)]
    breed_ebv = breed_ebv.assign(rel=(breed_ebv[e]/100)**2)
    breed_ebv = breed_ebv.dropna(subset=[d,'rel'])
    breed_ebv = breed_ebv[(breed_ebv.rel >= g)]
    print("Organizing the dEBV data...")
    data_debv = data_debv.rename({'ID': 'lpid'}, axis=1)
    data_ebv_debv = pd.merge(breed_ebv, data_debv, 
                       on='lpid', 
                       how='inner')
    data_ebv_debv = data_ebv_debv[['lpid', 'breed', 'flock', d,'rel','dEBV']]
    data_ebv_debv = data_ebv_debv[~data_ebv_debv.isin([np.nan, np.inf, -np.inf]).any(1)]
    print("Writing the output...")
    data_ebv_debv.to_csv(h, index=None, header=None, sep=" ", mode="w")
    print("END!!!")
    print("Good Luck!")


#------------------------------------------------------------------------------------------------------#

#-------------------------------- EXECUTION -----------------------------------------------------------#
# setting the directory
# if in purdue use below
#path = "C:/Users/araujoa/OneDrive/Documentos/Arquivos/Doutorado UESB/work_purdue_1/PhD_Thesis/NSIP_Chapter_4/data_sets/ped_perf"
# if in home use below
path = "C:/Users/Usuario/Desktop/André/OneDrive/Documentos/Arquivos/Doutorado UESB/work_purdue_1/PhD_Thesis/NSIP_Chapter_4/data_sets/ped_perf"
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
# Parameter e: breed code
# Parameter f: minimum reliability
# Parameter g: name of the output file  
prep_dEBV_data('Ramb EBV accuracies.xlsx','EBV & accuracies','dEBVs_bwt.csv','bwt','bwt_acc',612,0.0625,'dEBVs_bwt.txt')

print('                  ')
print('Trait 2 running...')
print('                  ')
# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: breed code
# Parameter f: minimum reliability
# Parameter g: name of the output file  
prep_dEBV_data('Ramb EBV accuracies.xlsx','EBV & accuracies','dEBVs_pwt.csv','pwt','pwt_acc',612,0.0625,'dEBVs_pwt.txt')

print('                  ')
print('Trait 3 running...')
print('                  ')
# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: breed code
# Parameter f: minimum reliability
# Parameter g: name of the output file  
prep_dEBV_data('Ramb EBV accuracies.xlsx','EBV & accuracies','dEBVs_ywt.csv','ywt','ywt_acc',612,0.0625,'dEBVs_ywt.txt')

print('                  ')
print('Trait 4 running...')
print('                  ')
# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: breed code
# Parameter f: minimum reliability
# Parameter g: name of the output file  
prep_dEBV_data('Ramb EBV accuracies.xlsx','EBV & accuracies','dEBVs_ygfw.csv','ygfw','ygfw_acc',612,0.0625,'dEBVs_ygfw.txt')

print('                  ')
print('Trait 5 running...')
print('                  ')
# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: breed code
# Parameter f: minimum reliability
# Parameter g: name of the output file  
prep_dEBV_data('Ramb EBV accuracies.xlsx','EBV & accuracies','dEBVs_yfd.csv','yfd','yfd_acc',612,0.0625,'dEBVs_yfd.txt')

print('                  ')
print('Trait 6 running...')
print('                  ')
# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: breed code
# Parameter f: minimum reliability
# Parameter g: name of the output file  
prep_dEBV_data('Ramb EBV accuracies.xlsx','EBV & accuracies','dEBVs_nlb.csv','nlb','nlb_acc',612,0.0625,'dEBVs_nlb.txt')

print('                  ')
print('Trait 7 running...')
print('                  ')
# Function to prepare the NSIP performance data from excel sheets
# Parameter a: name of the excel file ebv data
# Parameter b: name of the excel sheet ebv data
# Parameter c: name of trait_ebv
# Parameter d: name of trait_acc
# Parameter e: breed code
# Parameter f: minimum reliability
# Parameter g: name of the output file  
prep_dEBV_data('Ramb EBV accuracies.xlsx','EBV & accuracies','dEBVs_nlw.csv','nlw','nlw_acc',612,0.0625,'dEBVs_nlw.txt')


print('                  ')
print('END !!!')
print('Good Luck !!!')
print('finished ' + datetime.now().strftime("%m/%d/%Y %H:%M:%S"))
print('                  ')
####################################################################################################################