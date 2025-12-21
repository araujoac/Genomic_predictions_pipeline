# Genomic_predictions_pipeline
Code used to make genetic and genomic predictions with BLUPf90 softwares and validate using LR method used in Araujo et al. (2021; https://doi.org/10.1111%2Fjbg.12748)

If you would like to refence anything used from this code, please use the above reference. 

The pipeline assumes that the variances components are known, so make the estimates first.

The pipeline was originaly proposed to work with adjusted phenotypes in a single trait animal model, wich means that the 1'u vector is the fixed effect and the order of the addtive genetic effect after renumbering with renumf90 is 2. Other information for the inputs for renumf90 can be seen in the RENUM_FUNCTION. For other situations make the appropiate changes in the renum, calc_TA_EBV, and calc_TA_GEBV functions. All other traits analized in Araujo et al. (2021) were analized with same code except by the RENUM_FUNCTION, which was made appropriate for the trait model.
           
A grid search for best alpha and betha (blending of G and A matrices) and tau and omega (scaling of genomic a pedigree info in the H matrix) may also be done. In cases where the above grid rearches are not required you can just use the desired values.

You must check the functions preGS, blup, and blup_ssGBLUP to make sure that all the options are apropriate to your situation.

Changing the other formulas beyong the ones with options to BLUPf90 softwares may cause improper use. In this case contact the author. 

The pipeline takes a parameter file with the information needed to run, make sure you have all files in the specified format. Example of the paramater file:
  
 RENUM_FUNCTION  # function with the renumf90 parameter file to run   
 renum_add_mat  
 TRAIT # trait under evaluation, it is going to be informed in the RENUM_FUNCTION and added to some of the outputs  
 bwt  
 ADDITIVE_VARIANCE   # addtive genetic variance previously estimated  
 0.0854  
 PHENOTYPE_FILE # file with the phenotypes to run the renumf90 program in the RENUM_FUNCTION  
 bwt_adj.txt  
 PEDIGREE_FILE   # file with the pedigree to run the renumf90 program in the RENUM_FUNCTION  
 rambuoilet_ped.txt  
 TRUNCATION_DATE    # Truncation date (yyyymmdd) to devide the data sets in whole and partial (LR method)  
 20160421
 PED_FLOCK_DOB_SEX_FILE # pedigree file with id, sire id, dam id, flock id, DOB (yyyy-mm-dd format), and sex columns  
 ped_flock_dob.txt  
 DOWNLOAD_PROGRAMS  # binary to inform the pipeline to download (1) or not (0) the BLUPf90 softwares  
 0
  
You should haveall files, the Genomic_predictions_pipeline.sh, and Genomic_predictions_pipeline_par_file.par
to run the Genomic_predictions_pipeline.sh pipeline. You can you the following command to run:
  echo Genomic_predictions_pipeline_par_file.par | ./Genomic_predictions_pipeline.sh > Genomic_predictions_pipeline.out

This code is free of use under the MIT Licence. 

For problems, bugs or doubts contact Andre C. Araujo by email: araujoa@purdue.edu | andre.araujo@acufastswine.com
