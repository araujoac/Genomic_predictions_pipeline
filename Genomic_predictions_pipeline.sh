#!/usr/bin/env bash
#
# Name: Genomic_predictions_pipeline.sh
# Function: Make genetic and genomic predictions with BLUPf90 softwares and validate using LR method used in Araujo et al. (2021; https://doi.org/10.1111%2Fjbg.12748). 
# Site: See this after
# Author: André C. Araujo
# Manteiner: André C. Araujo
#
#--------------------------------------------------------------------------------------------------#
# Details:  This pipeline is used to make genetic and genomic predictions using the BLUPf90 softwares, which was used in Araujo et al. (2021; https://doi.org/10.1111%2Fjbg.12748),
#           If you would like to refence anything used from this code, please use the above reference. 
#           The pipeline assumes that the variances components are known, so make this estimates first
#           The pipeline was originaly proposed to work with adjusted phenotypes in a single trait
#           animal model, wich means that the 1'u vector is the fixed effect and the order of the addtive 
#           genetic effect after renumbering with renumf90 is 2. Other information for the inputs for renumf90 can be seen in the RENUM_FUNCTION. For other situations make the appropiate
#           changes in the renum, calc_TA_EBV, and calc_TA_GEBV functions.
#           All other traits analized in Araujo et al. (2021) were analized with same code except by the RENUM_FUNCTION, which was made appropriate for the trait model.
#           A grid search for best alpha and betha (blending of G and A matrices) and tau and omega (scaling of genomic a pedigree info in the H matrix)
#           may also be done. In cases where the above grid rearches are not required you can just use the desired values.
#           You must check the functions preGS, blup, and blup_ssGBLUP to make sure that all the options are apropriate to your situation
#           Changing the other formulas beyong the ones with options to BLUPf90 softwares may cause improper use. In this case contact the author. 
# 
#           The pipeline takes a parameter file with the information needed to run, make sure you have
#           all files in the specified format. Example of the paramater file:
#           
#           RENUM_FUNCTION  # function with the renumf90 parameter file to run 
#           renum_add_mat
#           TRAIT # trait under evaluation, it is going to be informed in the RENUM_FUNCTION and added to some of the outputs
#           bwt
#           ADDITIVE_VARIANCE   # addtive genetic variance previously estimated
#           0.0854
#           PHENOTYPE_FILE # file with the phenotypes to run the renumf90 program in the RENUM_FUNCTION
#           bwt_adj.txt
#           PEDIGREE_FILE   # file with the pedigree to run the renumf90 program in the RENUM_FUNCTION
#           rambuilet_ped.txt
#           TRUNCATION_DATE    # Truncation date (yyyymmdd) to devide the data sets in whole and partial (LR method)
#           20160421
#           PED_FLOCK_DOB_SEX_FILE # pedigree file with id, sire id, dam id, flock id, DOB (yyyy-mm-dd format), and sex columns 
#           ped_flock_dob.txt
#           DOWNLOAD_PROGRAMS  # binary to inform the pipeline to download (1) or not (0) the BLUPf90 softwares
#           0
# 
#           You should haveall files, the Genomic_predictions_pipeline.sh, and Genomic_predictions_pipeline_par_file.par
#           to run the Genomic_predictions_pipeline.sh pipeline. You can you the following command to run:
#            
#           echo Genomic_predictions_pipeline_par_file.par | ./Genomic_predictions_pipeline.sh > Genomic_predictions_pipeline.out
#
#           This code is free of use under the MIT Licence. 
#
#           FFor problems, bugs or doubts contact Andre C. Araujo by email: araujoa@purdue.edu | andre.araujo@acufastswine.com
#
#
#--------------------------------------------------------------------------------------------------#
# Change log:
# Criation: v1.0 in 12/01/2021
# Editions: 
#--------------------------------------------------------------------------------------------------#
# Tested in: bash 4.4.19(1)
#--------------------------------------------------------------------------------------------------#
# Acknowledges: Brito's Lab - Purdue University, Ron Lewis - UNL, BLUFf90 developers
#--------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------#

#------------------------------ VARIABLES ---------------------------------------------------------#
VERSION="v1.0"
# variable that controls the downlad_program function, 0 is no download and 1 means download
DOWNLOAD_PROGRAMS=0
#--------------------------------------------------------------------------------------------------#

#-------------------------------- FUNCTIONS -----------------------------------------------------------#

# Function with the presentation of the program
present() {
echo "******************************************************************
** Genomic_predictions_pipeline_NSIP.sh: $VERSION               
**                                                              
**             Beginning of run !!!                             
**                                                              
** Author: André C. Araujo                                      
** Contact: araujoa@purdue                                                                  
**                                                                                    
**                                                              
**              Good Luck  !!!                                  
**                                                              
**  Started at: "$(date)"                                       
******************************************************************
"
}


# Function to download the aplication programs (step 1 EXECUTION)
download_programs() {
echo "***************************************"
echo "* Downloading programs ..."
echo "***************************************"
# Links to download the BLUPF90 programs, the chmod command is to active the program in the folder
wget http://nce.ads.uga.edu/html/projects/programs/Linux/64bit/renumf90
chmod 777 renumf90
wget http://nce.ads.uga.edu/html/projects/programs/Linux/64bit/airemlf90
chmod 777 airemlf90
wget http://nce.ads.uga.edu/html/projects/programs/Linux/64bit/gibbs2f90
chmod 777 gibbs2f90
wget http://nce.ads.uga.edu/html/projects/programs/Linux/64bit/gibbs3f90
chmod 777 gibbs3f90
wget http://nce.ads.uga.edu/html/projects/programs/Linux/64bit/postgibbsf90
chmod 777 postgibbsf90
wget http://nce.ads.uga.edu/html/projects/programs/Linux/64bit/preGSf90
chmod 777 preGSf90
wget http://nce.ads.uga.edu/html/projects/programs/Linux/64bit/blupf90
chmod 777 blupf90
wget http://nce.ads.uga.edu/html/projects/programs/Linux/64bit/postGSf90
chmod 777 postGSf90
wget http://nce.ads.uga.edu/html/projects/programs/Linux/64bit/predictf90
chmod 777 predictf90
}


#
# Renum function for bwt - Birth Weight
#
# Function to create the parameter file and run renumf90 with maternal effects
# Parameter $1: phenotype data
# Model: Y = 1u + Za + Wm + Cc + e
renum_add_mat_bwt() {
# inserting the vector to account the mean
awk '{print $1,$2,$3,$4,$5,$6,$7=1}' $1 > $1.temp
rm $1
mv $1.temp $1
data_file=$1
####################################################
##
## Create the parameterfile to renumf90
##
####################################################

echo "DATAFILE
$data_file
TRAITS
6
FIELDS_PASSED TO OUTPUT

WEIGHT(S)

RESIDUAL_VARIANCE
0.3721
EFFECT
7 cross alpha #fixed overal mean
EFFECT
1 cross alpha #animal
RANDOM
animal
OPTIONAL
mat mpe
FILE
rambuilet_ped.txt
FILE_POS
1 2 3 0 0
SNP_FILE
snps.txt
PED_DEPTH
0
INBREEDING
pedigree
(CO)VARIANCES
0.0854 0
0 0.0915
(CO)VARIANCES_MPE
0.061
OPTION missing -999" > renum.par

####################################################
##
## Runing renumf90
##
####################################################

echo renum.par | ./renumf90 | tee renumf90.log
}


#
# BLUP function to predict EBV
#
# Function to predict breeding values with A matrix by removing SNP information
blup() {
#######################################################################################
# 
#     Solving MME using BLUP                        
#
#######################################################################################

# Remove the information related to markers for estimation of variance components and create the parameter file for airemlf90
sed -e '/SNP_FILE/d' renf90.par |
sed -e '/snps.txt/d'  > blup.par

# Include the options in renf90.par for blupf90
echo "OPTION conv_crit 1d-12" >> blup.par
echo "OPTION maxrounds 1000" >> blup.par
echo "OPTION sol se" >> blup.par
echo "OPTION residual" >> blup.par
echo "OPTION use_yams" >> blup.par

# Run airemlf90
echo blup.par | ./blupf90 | tee blup.log
cp solutions solutions_EBV.txt
cp yhat_residual yhat_residual_EBV.txt
}


#
# PREGS function to make quality control and save desired files
#
# Function to create the matrices for genomic evaluatiuon
# Parameter $1= tau value
# Parameter $2= omega value
# Parameter $3= alpha value
# Parameter $4= beta value
preGS() {
tau=$1
omega=$2
alpha=$3
beta=$4
#######################################################################################
# 
#     Create the matrices for genomic evaluation                       
#
######################################################################################
sed -i '/OPTION missing -999/d' renf90.par
# including the desired options to run preGSf90
echo "OPTION saveCleanSNPs" >> renf90.par
echo "OPTION map_file map.txt" >> renf90.par
echo "OPTION verify_parentage 1" >> renf90.par
echo "OPTION saveGInverse" >> renf90.par
echo "OPTION saveA22Inverse" >> renf90.par
echo "OPTION saveGOrig" >> renf90.par
echo "OPTION saveHinv" >> renf90.par
echo "OPTION saveAinv" >> renf90.par
echo "OPTION saveHinvOrig" >> renf90.par
echo "OPTION saveAinvOrig" >> renf90.par
echo "OPTION TauOmega $tau $omega" >> renf90.par
echo "OPTION AlphaBeta $alpha $beta" >> renf90.par

# Running preGSf90
echo renf90.par | ./preGSf90 | tee preGSf90.log
}


#
# BLUP function to predict GEBV
#
# Function to predict breeding values using ssGBLUP 
blup_ssGBLUP() {
#######################################################################################
# 
#     Solving MME using BLUP under ssGBLUP                       
#
######################################################################################
sed -e '/OPTION*/d' renf90.par > blup_ssGWAS.par

# Include the options in blup.par for blupf90
# replace the add_animal covariance type by user_file
sed -i "s/add_an_upginb/user_file/g" blup_ssGWAS.par
# replace the ped files by user_file relationship matrices, Hvinv matrices
sed -i "s/renadd02.ped/Hinv.txt/g" blup_ssGWAS.par

# Include the options in renf90.par for blupf90
echo "OPTION SNP_file snps.txt_clean" >> blup_ssGWAS.par
echo "OPTION no_quality_control" >> blup_ssGWAS.par
echo "OPTION conv_crit 1d-12" >> blup_ssGWAS.par
echo "OPTION maxrounds 1000" >> blup_ssGWAS.par
echo "OPTION use_yams" >> blup_ssGWAS.par
echo "OPTION sol se" >> blup_ssGWAS.par
echo "OPTION residual" >> blup_ssGWAS.par
echo "OPTION snp_p_value" >> blup_ssGWAS.par
echo "OPTION missing -999" >> blup_ssGWAS.par

# Run airemlf90
echo blup_ssGWAS.par | ./blupf90 | tee blup_ssGWAS.log
cp solutions solutions_GEBV.txt
cp yhat_residual yhat_residual_GEBV.txt
}


#
# POSTGS function 
#
# Function to run post genomic analyisis
# Parameter $1: Windows size to calculate the variances
postGS() {
win_size=$1
#######################################################################################
# 
#     post genomic analysis                       
#
######################################################################################
sed -e '/OPTION*/d' renf90.par > postgs.par

# Include the options in renf90.par for blupf90
echo "OPTION SNP_file snps.txt" >> postgs.par
echo "OPTION map_file map.txt" >> postgs.par
echo "OPTION no_quality_control" >> postgs.par
echo "OPTION snp_p_value" >> postgs.par
echo "OPTION readGInverse" >> postgs.par
echo "OPTION readA22Inverse" >> postgs.par
echo "OPTION which_weight 4" >> postgs.par
echo "OPTION windows_variance $win_size" >> postgs.par

# Running preGSf90
echo postgs.par | ./postGSf90 | tee postGSf90.log
}


# Function to calculate the theoretical accuracy form the blupf90 outputs for the EBVs
# Parameter $1: addtitive variance
# Parameter $2: number of the additive genetic effect from the solutions file
# Parameter $3: name of the output file
# Input files: solutions with se column and renf90.inb must be in the folder
# Output file: a file with five columns and header
calc_TA_EBV() {
add_var=$1
add_eff=$2
awk -v add_eff="$add_eff" '$2 == add_eff {print}' solutions > add_eff.temp # filter the additive effects from solutions
join -1 3 -2 3 <(sort -k3,3 add_eff.temp) <(sort -k3,3 renf90.inb) |       # merge the additive effects and A inb
awk -v add_var="$add_var" ' BEGIN {print "ID BV se inb TA"}
# calculate the theoretical accuracy
{$8=sqrt(1-($5**2/((1+$7)*add_var)));
if($8 ~/nan/) {print $6,$4,$5,$7,$8=0}
else {print $6,$4,$5,$7,$8} 
}' > $3 2> /dev/null
rm add_eff.temp
}


# Function to calculate the theoretical accuracy (TA) form the blupf90 outputs for the GEBVs
# Parameter $1: addtitive variance
# Parameter $2: number of the additive genetic effect from the solutions file
# Parameter $3: name of the output file
# Input files: solutions with se column, G_Orig.txt and renf90.inb must be in the folder
# Output file: on file with six columns and header having the TA, three inb files
calc_TA_GEBV() {
add_var=$1
add_eff=$2
awk -v add_eff="$add_eff" '$2 == add_eff {print}' solutions > add_eff.temp # finter the additive effects from solutions
awk '$1 == $2 {$4=$3-1; print}' G_Orig.txt > renf90_G_.inb.temp            # Calculate the G inb
join -1 1 -2 1 <(sort -k1,1 renf90.inb) <(sort -k1,1 renf90_G_.inb.temp) | # merge the full A inb with G inb 
awk '{print $1,$NF,$3,$4="Genotyped"}' > renf90_G_.inb                     # extract the file with G inb and solution codes
join -1 1 -2 1 <(sort -k1,1 renf90_G_.inb) <(sort -k1,1 renf90.inb) -v 2 | # merge keeping only the ungenotyped ids
awk '{$4="Non-genotyped"; print}' > renf90_A_non_G_.inb
cat renf90_A_non_G_.inb renf90_G_.inb > renf90_H_.inb                      # create the H inb file
join -1 3 -2 3 <(sort -k3,3 add_eff.temp) <(sort -k3,3 renf90_H_.inb) |    # merge the additive effects and H inb
awk -v add_var="$add_var" 'BEGIN {print "ID BV se inb TA Genotypes"}
# calculate the theoretical accuracy
{$9=sqrt(1-($5**2/((1+$7)*add_var))); 
if($9 ~/nan/) {print $6,$4,$5,$7,$9=0,$8}
else {print $6,$4,$5,$7,$9,$8}
}' > $3 2> /dev/null
rm  add_eff.temp renf90_G_.inb.temp
}


# Function to cretate a file with genotyped individuals with phenotypes or sires or dams not phenotyped, but with phenotyped progeny
# Parameter $1: name of a pedigree file with id, sire id, dam id, flock id, DOB (yyyy-mm-dd format), and sex columns
# Parameter $2: name snp file with id and snp string column (or only the id column)
# Parameter $3: name of the phenotype file with the id in the first column (other columns may be present, but not used)
# Output 1: ped_flock_dob_geno.txt, with the pedigree for genotyped ids including flock and DOB
# Output 2: ped_flock_dob_pheno.txt, with the pedigree for phenotyped ids including flock and DOB
# Output 3: ped_flock_dob_geno_pheno_or_p_pheno.txt, with genotyped ids that matches the function purposes including a column to classify ids
create_genotypes_to_divide_training_and_validation() {
echo "*** create_genotypes_to_divide_training_and_validation function running...
*** Author: Andre C. Araujo
*** Contact: araujoa@purdue.edu
*** Started: "$(date)""
    # cerat ped_flock_dob file for genotyped animals
    join -1 1 -2 1 <(awk '{$6=""; print}' $1 | sort -k1,1) <(awk '{print $1}' $2 | sort -k1,1) > ped_flock_dob_geno.txt
    # create ped_flock_dob file for phenotyped animals
    join -1 1 -2 1 <(awk '{$6=""; print}' $1 | sort -k1,1) <(awk '{print $1}' $3 | sort -k1,1) | sort -k1,1 -u > ped_flock_dob_pheno.txt
    # create the ped_flock_dob file for genotyped individuals that have phenotypes
    awk '
      FNR==NR{
        i[$1];j[$2];k[$3];next
      }
      ($1 in i) { 
        $6="phen"; print 
      }
    ' ped_flock_dob_pheno.txt ped_flock_dob_geno.txt > ped_flock_dob_geno_pheno.temp
    # create the ped_flock_dob file for sires with no phenotypes but with progeny with phenotypes
    awk '
      FNR==NR{
        i[$1];j[$2];k[$3];next
      }
      ($1 in j) && !($1 in i) { 
        $6="sire_no_phen"; print 
      }
    ' ped_flock_dob_pheno.txt ped_flock_dob_geno.txt > ped_flock_dob_sire_geno.temp
    # create the ped_flock_dob file for dams with no phenotypes but with progeny with phenotypes
    awk '
      FNR==NR{
        i[$1];j[$2];k[$3];next
      }
      ($1 in k) && !($1 in i) { 
        $6="dam_no_phen"; print 
      }
    ' ped_flock_dob_pheno.txt ped_flock_dob_geno.txt > ped_flock_dob_dam_geno.temp
    # statistics idviduals 
    phen=$(awk '{print $1}' ped_flock_dob_geno_pheno.temp | wc -l)
    sire_no_phen=$(awk '{print $1}' ped_flock_dob_sire_geno.temp | wc -l)
    dam_no_phen=$(awk '{print $1}' ped_flock_dob_dam_geno.temp | wc -l)
    ids_ped=$(awk '{print $1}' $1 | wc -l)
    ids_genotyped=$(awk '{print $1}' $2 | wc -l)
    ids_phen=$(awk '{print $1}' $3 | sort -k1,1 -u | wc -l)
    ids_genotyped_in_ped=$(awk '{print $1}' ped_flock_dob_geno.txt | wc -l)
    ids_phenotyped_in_ped=$(awk '{print $1}' ped_flock_dob_pheno.txt | wc -l)
# print statistics
echo "
Number of individuals with parent information: $ids_ped

Genotyped inividuals in the SNP file: $ids_genotyped
Genotyped inividuals with parent information: $ids_genotyped_in_ped
Genotyped inividuals with phenotype and parent information (phen): $phen
Genotyped sires with no phenotypes (sire_no_phen): $sire_no_phen
Genotyped dams with no phenotypes (dam_no_phen): $dam_no_phen
Total genotyped individuals cosidered or further analysis (phen+sire_no_phen+dam_no_phen): $(( $phen+$sire_no_phen+$dam_no_phen ))

Number of individuals with phenotype: $ids_phen
Number of inividuals with phenotype and parent information: $ids_phenotyped_in_ped"
    # create the output
    cat ped_flock_dob_geno_pheno.temp ped_flock_dob_sire_geno.temp ped_flock_dob_dam_geno.temp > ped_flock_dob_geno_pheno_or_p_pheno.txt
    rm ped_flock_dob_geno_pheno.temp ped_flock_dob_sire_geno.temp ped_flock_dob_dam_geno.temp
echo "
Files ped_flock_dob_geno.txt, ped_flock_dob_pheno.txt, and ped_flock_dob_geno_pheno_or_p_pheno.txt written.

*** Finished: "$(date)""
}


# Function to create a file with calling training and validation individuals
# No header in the input 
# Parameter $1: truncation date in the yyyymmdd WITH NO SEPARATORS
# Parameter $2: name of the input file with extension (Output 3 from create_genotypes_to_divide_training_and_validation function)
# Parameter $3: name of the theoretical accuracy (TA) file with extension (from calc_TA_EBV function)
# Parameter $4: name of the phenotype file with extension 
# Parameter $5: name of the input file with extension (Output 2 from create_genotypes_to_divide_training_and_validation function)
# Output 1: file named training_validation.txt similiar to input, adding columns specifying trainig or validation for each id and TA   
# Output 2: file named partial_data.txt similiar to the phenotype file, removing all phenotypes recorded after the truncation date
create_training_validation() {
echo "*** create_training_validation function running...
*** Author: Andre C. Araujo
*** Contact: araujoa@purdue.edu
*** Started: "$(date)""
    yyyymmdd=$1
    # create training and validation
    awk '{gsub("-","",$5); print}' $2 | 
    sort -k5,5 -nr |
    awk -v yyyymmdd="$yyyymmdd" '{if($5 < yyyymmdd) {$7="training"; print} else {$7="validation"; print}}' |
    sort -k5,5 -nr > training_validation.temp
    # include the theoretical accuracy in the column $8
    join -1 1 -2 1 <(sort -k1,1 training_validation.temp) <(sort -k1,1 $3) | 
    awk '{print $1,$2,$3,$4,$5,$6,$7,$NF}' > training_validation.txt
    rm training_validation.temp
    # creating the partial dataset removing phenotypes recorded after yyyymmdd
    join -1 1 -2 1 <(sort -k1,1 $4) <(sort -k1,1 $5) |
    awk -v yyyymmdd="$yyyymmdd" '{gsub("-","",$NF); if($NF < yyyymmdd) print $1,$2,$3,$4,$5,$6,$7}' > partial_data.txt
    # statistics about training and validation
    training=$(awk '$7 == "training" {print $1}' training_validation.txt | wc -l)
    validation=$(awk '$7 == "validation" {print $1}' training_validation.txt | wc -l)
    total=$(awk '{print $1}' $2 | wc -l)
    percent_validation=$(awk -v training="$training" -v validation="$validation" '
    BEGIN {total=training+validation; percent_validation=(validation/total)*100; print percent_validation}')
    TA_training=$(awk -v training="$training" '$7 == "training" {sum+=$8} 
    END {print sum/training}' training_validation.txt)
    TA_validation=$(awk -v validation="$validation" '$7 == "validation" {sum+=$8} 
    END {print sum/validation}' training_validation.txt)
    Overall_TA=$(awk '{sum+=$8} END {print sum/NR}' training_validation.txt)
    min_TA_training=$(awk '$7 == "training" {print $8}' training_validation.txt | sort -n | head -1)
    max_TA_training=$(awk '$7 == "training" {print $8}' training_validation.txt | sort -nr | head -1)
    min_TA_validation=$(awk '$7 == "validation" {print $8}' training_validation.txt | sort -n | head -1)
    max_TA_validation=$(awk '$7 == "validation" {print $8}' training_validation.txt | sort -nr | head -1)
    min_Overall_TA_training=$(awk '{print $8}' training_validation.txt | sort -n | head -1)
    max_Overall_TA_training=$(awk '{print $8}' training_validation.txt | sort -nr | head -1)
    pheno_whole_dataset=$(awk '{print $1}' $4 | wc -l)
    pheno_partial_dataset=$(awk '{print $1}' partial_data.txt | wc -l)
# print statistics
echo "
Total individuals: $total
Training individuals: $training
Validation individuals: $validation
Percent validation individuals: $percent_validation%

Number of records whole dataset: $pheno_whole_dataset
Number of records partial dataset: $pheno_partial_dataset

Average (min-max) Theoretical Accuracy for training individuals: $TA_training ($min_TA_training-$max_TA_training)
Average (min-max) Theoretical Accuracy for validation individuals: $TA_validation ($min_TA_validation-$max_TA_validation)
Average (min-max) Theoretical Accuracy for all individuals: $Overall_TA ($min_Overall_TA_training-$max_Overall_TA_training)

Files training_validation.txt and partial_data.txt writen.
"
echo "*** Finished: "$(date)""
}


# Function run a grid research for the tau and omega parameters in the H inverse creation
# Parameter $1: alpha parameter
# Parameter $2: tau grid
# Parameter $3: omega grid  
tau_omega_grid_research() {
# create a folder to run the different tau and omega values for the alpha equal to 0.95
alpha=$1
rm -rf diff_t_o_$(echo $alpha)_G
mkdir diff_t_o_$(echo $alpha)_G
cp *f90 "$phenotype" EBV_TA_$(echo $trait).txt $(echo $pedigree) snps.txt map.txt diff_t_o_$(echo $alpha)_G/
cd diff_t_o_$(echo $alpha)_G
beta=$(awk -v alpha="$alpha" 'BEGIN {omega=1-alpha; print omega}' )
# Running the genomic predictions over different tau and omega
for i in $2; do  # tau values grid
    for j in $3; do  # omega values grid
        rm -rf genomic_$(echo $alpha)_G_t_$(echo $i)_o_$(echo $j)
        mkdir genomic_$(echo $alpha)_G_t_$(echo $i)_o_$(echo $j)
        cp *f90 "$phenotype" EBV_TA_$(echo $trait).txt $(echo $pedigree) snps.txt map.txt genomic_$(echo $alpha)_G_t_$(echo $i)_o_$(echo $j)/
        cd genomic_$(echo $alpha)_G_t_$(echo $i)_o_$(echo $j)

        echo " "
        echo "Tau $i and Omega $j"
        pwd
        echo " "

        # Function to create the parameter file and run renumf90
        $(echo $renum)_$(echo $trait) "$phenotype"

        # Function to create the matrices for genomic evaluatiuon
        # Parameter $1= tau value
        # Parameter $2= omega value
        # Parameter $3= alpha value
        # Parameter $4= beta value
        preGS "$i" "$j" "$alpha" "$beta"

        # function to predict breeding values using ssGBLUP 
        blup_ssGBLUP

        # Function to calculate the theoretical accuracy (TA) form the blupf90 outputs for the GEBVs
        # Parameter $1: addtitive variance
        # Parameter $2: number of the additive genetic effect from the solutions file
        # Parameter $3: name of the output file
        # Input files: solutions with se column, G_Orig.txt and renf90.inb must be in the folder
        # Output file: on file with six columns and header having the TA, three inb files
        calc_TA_GEBV "$add_var" "2" "GEBV_TA_"$(echo $trait)".txt"

        # Function to merge the files with GEBV, EBV, Theoretical Accuracy, and Iinbreeding abd 
		# Parameter $1: trait
		# Parameter $2: trait heritability
		marge_gebv_ebv_inb_ta_files "$trait" "h"

        cd ../
    done
done
cd ../
}


# Function run a grid research for the tau and omega parameters in the H inverse creation (partial data set)
# Parameter $1: alpha parameter
# Parameter $2: tau grid
# Parameter $3: omega grid 
tau_omega_grid_research_partial_data() {
# create a folder to run the different tau and omega values for the alpha equal to 0.95
alpha=$1
rm -rf diff_t_o_$(echo $alpha)_G_partial_data
mkdir diff_t_o_$(echo $alpha)_G_partial_data
cp *f90 partial_data.txt $(echo $pedigree) snps.txt map.txt diff_t_o_$(echo $alpha)_G_partial_data/
cd diff_t_o_$(echo $alpha)_G_partial_data
beta=$(awk -v alpha="$alpha" 'BEGIN {omega=1-alpha; print omega}' )
# Running the genomic predictions over different tau and omega
for i in $2; do  # tau values grid
    for j in $3; do  # omega values grid
        rm -rf genomic_$(echo $alpha)_G_t_$(echo $i)_o_$(echo $j)
        mkdir genomic_$(echo $alpha)_G_t_$(echo $i)_o_$(echo $j)
        cp *f90 partial_data.txt $(echo $pedigree) snps.txt map.txt genomic_$(echo $alpha)_G_t_$(echo $i)_o_$(echo $j)/
        cd genomic_$(echo $alpha)_G_t_$(echo $i)_o_$(echo $j)

        echo " "
        echo "Tau $i and Omega $j"
        pwd
        echo " "

        # Function to create the parameter file and run renumf90
        # Parameter $1: phenotype data
        $(echo $renum)_$(echo $trait) "partial_data.txt"

        # Function to create the matrices for genomic evaluatiuon
        # Parameter $1= tau value
        # Parameter $2= omega value
        # Parameter $3= alpha value
        # Parameter $4= beta value
        preGS "$i" "$j" "$alpha" "$beta"

        # function to predict breeding values using ssGBLUP 
        blup_ssGBLUP

        # Function to calculate the theoretical accuracy (TA) form the blupf90 outputs for the GEBVs
        # Parameter $1: addtitive variance
        # Parameter $2: number of the additive genetic effect from the solutions file
        # Parameter $3: name of the output file
        # Input files: solutions with se column, G_Orig.txt and renf90.inb must be in the folder
        # Output file: on file with six columns and header having the TA, three inb files
        calc_TA_GEBV "$add_var" "2" "GEBV_TA_partial_data.txt"

        cd ../
    done
done
cd ../
}


# Function run a grid research for the tau and omega parameters in the H inverse creation (partial data set)
# Parameter $1: alpha parameter
# Parameter $2: tau grid
# Parameter $3: omega grid 
# Parameter $4: additive variance 
apply_LR_to_tau_omega_partial_data_resuls() {
# create a folder to run the different tau and omega values for the alpha equal to 0.95
alpha=$1
add_var=$4

echo " "
echo "apply_LR_to_tau_omega_partial_data_resuls function running alpha $alpha"

# going to the whole data set folder
cd diff_t_o_$(echo $alpha)_G
# going to the different folders with whole data results for different tau and omega
for i in $2; do  # tau values grid
    for j in $3; do  # omega values grid
        cd genomic_$(echo $alpha)_G_t_$(echo $i)_o_$(echo $j)
        echo " "
        echo "Tau $i and Omega $j"
        pwd
        cp GEBV_TA_$(echo $trait).txt GEBV_TA_$(echo $alpha)_G_t_$(echo $i)_o_$(echo $j)_whole_dataset.txt
        cp GEBV_TA_$(echo $alpha)_G_t_$(echo $i)_o_$(echo $j)_whole_dataset.txt ../../
        cd ../
    done
done
cd ../
# going to the partial data set folder
cd diff_t_o_$(echo $alpha)_G_partial_data
# going to the different folders with partial data results for different tau and omega
for i in $2; do  # tau values grid
    for j in $3; do  # omega values grid
        cd genomic_$(echo $alpha)_G_t_$(echo $i)_o_$(echo $j)
        echo " "
        echo "Tau $i and Omega $j"
        pwd
        cp GEBV_TA_partial_data.txt GEBV_TA_$(echo $alpha)_G_t_$(echo $i)_o_$(echo $j)_partial_dataset.txt
        cp GEBV_TA_$(echo $alpha)_G_t_$(echo $i)_o_$(echo $j)_partial_dataset.txt ../../
        cd ../
    done
done
cd ../

# create a function to calculate the LR statistics
rm -rf LR_results_t_o_grid_research_$(echo $alpha)_G 
mkdir LR_results_t_o_grid_research_$(echo $alpha)_G
cp training_validation.txt EBV_TA_$(echo $trait).txt EBV_TA_partial.txt lr.R LR_results_t_o_grid_research_$(echo $alpha)_G/
mv *_whole_dataset.txt LR_results_t_o_grid_research_$(echo $alpha)_G/
mv *_partial_dataset.txt LR_results_t_o_grid_research_$(echo $alpha)_G/
cd LR_results_t_o_grid_research_$(echo $alpha)_G/

# going to the different folders with partial data results for different tau and omega
for i in $2; do  # tau values grid
    for j in $3; do  # omega values grid
        echo " "
        echo "Getting LR results with alpha $alpha for Tau $i and Omega $j"
        # merge the GEBV whole with the training_and_validation filetering only the validation ids
        join -1 1 -2 1 <(awk '{print $1,$2,$3,$4,$5}' GEBV_TA_$(echo $alpha)_G_t_$(echo $i)_o_$(echo $j)_whole_dataset.txt | sort -k1,1) \
        <(awk '$7 == "validation" {print $1,$6}' training_validation.txt | sort -k1,1) > validation.temp
        # merge previous file with EBV whole 
        join -1 1 -2 1 <(awk '{print $1,$6,$4,$2,$3,$5}' validation.temp | sort -k1,1) \
        <(awk '{print $1,$2,$3,$4,$5}' EBV_TA_$(echo $trait).txt | sort -k1,1) |
        awk '{print $1,$2,$3,$9,$4,$5,$6,$7,$8,$10}' > validation_2.temp
        # merge previous file with GEBV partial
        join -1 1 -2 1 <(awk '{print}' validation_2.temp | sort -k1,1) \
        <(awk '{print $1,$2,$3,$5}' GEBV_TA_$(echo $alpha)_G_t_$(echo $i)_o_$(echo $j)_partial_dataset.txt | sort -k1,1) > validation_3.temp
        # merge previous file with EBV partial
        join -1 1 -2 1 <(awk '{print}' validation_3.temp | sort -k1,1) \
        <(awk '{print $1,$2,$3,$NF}' EBV_TA_partial.txt | sort -k1,1) |
        awk 'BEGIN {print "id phen_st G_inb A_inb GEBV_w GEBV_sd_w GEBV_TA_w EBV_w EBV_sd_w EBV_TA_w GEBV_p GEBV_sd_p GEBV_TA_p EBV_p EBV_sd_p EBV_TA_p"}
        {print}' > validation.txt
        rm validation.temp validation_2.temp validation_3.temp
        # applying the LR method
        sed -i "s/sigma2a/$add_var/g" lr.R
        Rscript lr.R
        rm validation.txt
        sed -i 's/"//g' lr_result.txt
        awk -v i="$i" -v j="$j" '{$23=i"_"j; print}' lr_result.txt > lr_result_$(echo $alpha)_G_t_$(echo $i)_o_$(echo $j).txt
        rm lr_result.txt
    done
done
echo "
LR results $trait:
"
cat lr_result_* |
awk 'BEGIN {print "t_o GEBVaccW GEBVbiasW GEBVdispW EBVaccW EBVbiasW EBVdispW GEBVTAW EBVTAW GEBVTAP EBVTAP GEBVseW EBVseW GEBVseP EBVseP"} 
{print $23,$1,$2,$3,$4,$5,$6,$7"("$8")",$9"("$10")",$11"("$12")",$13"("$14")",$15"("$16")",$17"("$18")",$19"("$20")",$21"("$22")"}' > lr_results.txt
cat lr_results.txt
echo " "
echo "END"
cd ../
}

# Function to merge the files with GEBV, EBV, Theoretical Accuracy, and Iinbreeding abd 
# Parameter $1: trait
# Parameter $2: trait heritability
marge_gebv_ebv_inb_ta_files(){
h=$2
join -1 1 -2 1 <(awk 'FNR > 1 {print}' GEBV_TA_$1.txt | sort -k1,1) \
<(awk 'FNR > 1 {print}' EBV_TA_$1.txt | sort -k1,1) |
sort -k1,1 -k7,7 -n |
awk -v h="$h" 'BEGIN {print "ID GEBV GEBV_se H_inb TA_GEBV ENP_GEBV Genotypes EBV EBV_se A_inb TA_EBV ENP_EBV"}
{   # calculate the Effective number of progeny for the GEBV and EBV ($11 and $12, respectively)
	$11=(($5**2)*((4-(h**2))/(h**2)))/(1-($5**2));
	$12=(($10**2)*((4-(h**2))/(h**2)))/(1-($10**2)); 
	print $1,$2,$3,$4,$5,$11,$6,$7,$8,$9,$10,$12
}' > gebv_ebv_inb_ta_$1.txt	
}

#------------------------------------------------------------------------------------------------------#

#-------------------------------- R CODES -------------------------------------------------------------#

# R code to apply the LR method
echo 'data=read.table("validation.txt", header=T)

add_var=sigma2a

# calculate the GEBV accuracy
ave_G_inb=mean(data$G_inb)
cov_GEBV_wp=cov(data$GEBV_w,data$GEBV_p)
GEBV_acc=sqrt(cov_GEBV_wp/((1-ave_G_inb)*add_var))
# calulate GEBV bias
GEBV_bias=mean(data$GEBV_p)-mean(data$GEBV_w)
# calculate the GEBV dispersion
GEBV_disp=cov_GEBV_wp/var(data$GEBV_p)

# calculate the EBV accuracy
ave_A_inb=mean(data$A_inb)
cov_EBV_wp=cov(data$EBV_w,data$EBV_p)
EBV_acc=sqrt(cov_EBV_wp/((1-ave_A_inb)*add_var))
# calulate GEBV bias
EBV_bias=mean(data$EBV_p)-mean(data$EBV_w)
# calculate the GEBV dispersion
EBV_disp=cov_EBV_wp/var(data$EBV_p)

# average Thoretical Accurracies
ave_GEBV_TA_w=mean(data$GEBV_TA_w)
ave_EBV_TA_w=mean(data$EBV_TA_w)
ave_GEBV_TA_p=mean(data$GEBV_TA_p)
ave_EBV_TA_p=mean(data$EBV_TA_p)

# average standard errors Accurracies
ave_GEBV_se_w=mean(data$GEBV_sd_w)
ave_EBV_se_w=mean(data$EBV_sd_w)
ave_GEBV_se_p=mean(data$GEBV_sd_p)
ave_EBV_se_p=mean(data$EBV_sd_p)

# standard deviations Thoretical Accurracies
sd_GEBV_TA_w=sd(data$GEBV_TA_w)
sd_EBV_TA_w=sd(data$EBV_TA_w)
sd_GEBV_TA_p=sd(data$GEBV_TA_p)
sd_EBV_TA_p=sd(data$EBV_TA_p)

# standard deviations standard errors Accurracies
sd_GEBV_se_w=sd(data$GEBV_sd_w)
sd_EBV_se_w=sd(data$EBV_sd_w)
sd_GEBV_se_p=sd(data$GEBV_sd_p)
sd_EBV_se_p=sd(data$EBV_sd_p)

lr_result=cbind(GEBV_acc,GEBV_bias,GEBV_disp,EBV_acc,EBV_bias,EBV_disp,
    ave_GEBV_TA_w,sd_GEBV_TA_w,ave_EBV_TA_w,sd_EBV_TA_w,
    ave_GEBV_TA_p,sd_GEBV_TA_p,ave_EBV_TA_p,sd_EBV_TA_p,
    ave_GEBV_se_w,sd_GEBV_se_w,ave_EBV_se_w,sd_EBV_se_w,
    ave_GEBV_se_p,sd_GEBV_se_p,ave_EBV_se_p,sd_EBV_se_p)
lr_result=format(round(lr_result, 3), nsmall = 3)
write.table(lr_result,"lr_result.txt", row.names=F, col.names=F)
q()
n' > lr.R


Rscript lr.R

#------------------------------------------------------------------------------------------------------#

#-------------------------------- EXECUTION -----------------------------------------------------------#

### 0) Presentation and settig the information to run the program

# 0.1) Function with the presentation of the program
present

# 0.2) Show basic information
echo 'Make sure the phenotype, pedigree, ped_flock_dob_sex, SNP, and map files in the folder before running.'
echo 'See prep_data_pipeline_NSIP.sh for more information.'
echo ' '
echo 'Make sure you have the appropriate renumf90 parameter files in the renum functions before running.'

echo '
Inform parfile:'
read parfile

# 0.3) read the input variables
renum=$(sed '2!d' $(echo $parfile)) # Extracting renum function
trait=$(sed '4!d' $(echo $parfile)) # Extracting trait
add_var=$(sed '6!d' $(echo $parfile)) # Extracting add_var
phenotype=$(sed '8!d' $(echo $parfile)) # Extracting phenotype
pedigree=$(sed '10!d' $(echo $parfile)) # Extracting pedigree
trunc_date=$(sed '12!d' $(echo $parfile)) # Extracting snp_file
ped_flock_dob_sex=$(sed '14!d' $(echo $parfile)) # Extracting ped_flock_dob_sex
DOWNLOAD_PROGRAMS=$(sed '16!d' $(echo $parfile)) # Extracting DOWNLOAD_PROGRAMS
echo ' '
echo "Renum function: $renum
Trait: $trait
Additive variance: $add_var
Phenotype: $phenotype
Pedigree: $pedigree
Truncation date: $trunc_date
ped_flock_dob_sex: $ped_flock_dob_sex
"


# 0.4) Download programns (optional step controled by DOWNLOAD_PROGRAMS variable, see VARIABLES section)
# Bash if to apply download_programs function
if [ $DOWNLOAD_PROGRAMS -eq 1 ]; then
    # Function to download the aplication programs
    download_programs
else
    echo 'BLUPf90 programs installed.'
fi

# 0.5) Set the modules
module load r
module load anaconda

#
## Analysis starts here
#
start=$SECONDS # take the time in the beginning (righ after read all functions)

echo '
*** 1) Analysis using pedigree
'
echo '** 1.1) Analysis using pedigree (A) matrix
'
echo '1.1.1) Whole data set
' 

# Run function to create the parameter file and run renumf90
$(echo $renum)_$(echo $trait) "$phenotype"

# Function to predict EBV with A matrix by removing SNP information
# Model: Y = 1u + Za + Wm + Cc + e
blup

# Function to calculate the theoretical accuracy form the blupf90 outputs
# Parameter $1: addtitive variance
# Parameter $2: number of the additive genetic effect from the solutions file
# Parameter $3: name of the pedigree file, with the year of birth at the last column
# Parameter $4: name of the output file
calc_TA_EBV "$add_var" "2" "EBV_TA_"$(echo $trait)".txt"

# Function to cretate a file with genotyped individuals with phenotypes or sires or dams not phenotyped, but with phenotyped progeny
# Parameter $1: name of a pedigree file with id, sire id, dam id, flock id, and DOB (yyyy-mm-dd format) columns
# Parameter $2: name snp file with id and snp string column (or only the id column)
# Parameter $3: name of the phenotype file with the id in the first column (other columns may be present, but not used)
# Output 1: ped_flock_dob_geno.txt, with the pedigree for genotyped ids including flock and DOB
# Output 2: ped_flock_dob_pheno.txt, with the pedigree for phenotyped ids including flock and DOB
# Output 3: ped_flock_dob_geno_pheno_or_p_pheno.txt, with genotyped ids that matches the function purposes including a column to classify ids
create_genotypes_to_divide_training_and_validation "$ped_flock_dob_sex" "snps.txt" "$phenotype"

# Function to create a file with calling training and validation individuals
# Input: file with 5 columns: id, sire id, dam id, flock, and date of birth (DOB) (yyyy-mm-dd format, no space or dash separated)
# No header in the input 
# Parameter $1: DOB in the yyyy-mm-dd WITH NO SEPARATORS
# Parameter $2: name of the input file with extension
# Output: file named training_validation.txt similiar to input, adding a column specifying trainig or validation for each id  
create_training_validation "$trunc_date" "ped_flock_dob_geno_pheno_or_p_pheno.txt" "EBV_TA_"$(echo $trait)".txt" \
"$phenotype" "ped_flock_dob_pheno.txt"


echo '1.1.1) Partial data set
' 

rm -rf A-BLUP_partial_data
mkdir A-BLUP_partial_data
cp *f90 partial_data.txt $(echo $pedigree) snps.txt map.txt A-BLUP_partial_data/
cd A-BLUP_partial_data/

# Function to create the parameter file and run renumf90
# Parameter $1: phenotype data
# Model: Y = 1u + Za + Wm + Cc + e
$(echo $renum)_$(echo $trait) "partial_data.txt"

# Run function to calculate EBV with A matrix by removing SNP information
blup

# Function to calculate the theoretical accuracy form the blupf90 outputs
# Parameter $1: addtitive variance
# Parameter $2: number of the additive genetic effect from the solutions file
# Parameter $3: name of the pedigree file, with the year of birth at the last column
# Parameter $4: name of the output file
calc_TA_EBV "$add_var" "2" "EBV_TA_partial.txt"

cp EBV_TA_partial.txt ../
cd ../

# 1.2) GCA analisis using pedigree (A) matrix
# Organizing the pedigree and pheno files to GCA pakage 
#join -1 10 -2 1 <(sort -k10,10 renadd02.ped) <(sort -k1,1 bwt_adj.txt) > pheno_ped.txt

#module load r
# Run the R code for GCA analysis
#R CMD BATCH GCA.R

###
echo '
*** 2) Analysis using genomic data
'
echo '** 2.1) Analysis using different alpha and beta and tau and omega parameters
'
echo '2.1.2) Analysis using whole data set
' 

# Function run a grid research for the tau and omega parameters in the H inverse creation
# Parameter $1: alpha parameter
# Parameter $2: tau grid
# Parameter $2: omega grid 
tau_omega_grid_research "0.95" "1.0" "1.0"


echo '
2.1.3) Analysis using partial data set
'

# Function run a grid research for the tau and omega parameters in the H inverse creation
# Parameter $1: alpha parameter
# Parameter $2: tau grid
# Parameter $2: omega grid 
tau_omega_grid_research_partial_data "0.95" "1.0" "1.0"


###
echo '
*** 3) Apply LR method to tau omega grid research
' 

# Function run a grid research for the tau and omega parameters in the H inverse creation (partial data set)
# Parameter $1: alpha parameter
# Parameter $2: tau grid
# Parameter $3: omega grid 
# Parameter $4: additive variance 
apply_LR_to_tau_omega_partial_data_resuls "0.95" "1.0 1.5 2.0" "0.4 0.7 1.0" "$add_var"
apply_LR_to_tau_omega_partial_data_resuls "0.50" "1.0 1.5 2.0" "0.4 0.7 1.0" "$add_var"


end=$SECONDS # take the time in the end
duration=$(( end - start )) # take the duration Genomic_predictions_duration_time

echo "Duration time in secons Genomic_predictions_pipeline_NSIP.sh : $duration " > Genomic_predictions_duration_time.txt

###
echo "
*** Finished: "$(date)"" 
  
#------------------------------------------------------------------------------------------------------#