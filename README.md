# SARS-CoV-2_90K_patient_study

## ABSTRACT 

### Purpose 

Glycoprotein 90K, encoded by the interferon-stimulated gene LGALS3BP, displays broad antiviral activity. It reduces HIV-1 infectivity by interfering with Env maturation and virion incorporation, and increases survival of Influenza A virus-infected mice via antiviral innate immune signaling. Its antiviral potential in SARS-CoV-2 infection remains largely unknown. 

### Methods 

Here, we analyzed the expression of 90K/LGALS3BP in 44 hospitalized COVID-19 patients at multiple levels. We quantified 90K protein concentrations in serum and PBMCs as well as LGALS3BP mRNA levels. Complementary, we analyzed two single cell RNA-sequencing datasets for expression of LGALS3BP in respiratory specimens and PBMCs from COVID-19 patients. Finally, we analyzed the potential of 90K to interfere with SARS-CoV-2 infection of HEK293T/ACE2, Calu-3 and Caco-2 cells using authentic virus. 

### Results 

90K protein serum concentrations were significantly elevated in COVID-19 patients compared to uninfected sex- and age-matched controls. Furthermore, PBMC-associated concentrations of 90K protein were overall reduced by SARS-CoV-2 infection in vivo, suggesting enhanced secretion into the extracellular space. Mining of published PBMC scRNA-seq datasets uncovered monocyte-specific induction of LGALS3BP mRNA expression in COVID-19 patients. In functional assays, neither 90K overexpression in susceptible cell lines nor exogenous addition of purified 90K consistently inhibited SARS-CoV-2 infection. 
Conclusion Our data suggests that 90K/LGALS3BP contributes to the global type I IFN response during SARS-CoV-2 infection in vivo without displaying detectable antiviral properties in vitro.

## OVERVIEW

This is the GitHub repository for the manuscript "90K/LGALS3BP Expression Is Upregulated in COVID-19 but may not Restrict SARS-CoV-2 Infection " by Bosquillon de Jarcy _et al._, 2023.  The code utilised for analysis of the patient data is deposited here. 

### Patient data analysis

All code required to perform analysis on patient data is provided in the folder 01_PatientDataAnalysis. The actual data is also available here, in the form of a .csv file (90K_data_all.csv). A brief overview of the column headings in this file and their meaning:

* pat_id: Patient ID; signifier used to uniquely identify the patients enrolled in the study.
* abnahme_id: Sampling ID; signifier used to uniquely identify a sampling timepoint for a patient.
* dpso: Days post symptom onset; the number of days that had passed since the onset of COVID-19 symptoms at the sampling timepoint.
* serum_90k: Concentration of 90k protein in the patient's serum, in μg/ml, measured by ELISA.
* healthy_serum_90k: Concentration of 90k protein in the matched healthy control's serum, when this is available for a patient, in μg/ml, measured by ELISA.
* protein_pbmc_90k: Concentration of 90k protein in PBMC cell lysates for the patientint, in μg/ml, measured by ELISA. "no" indicates that this was not measured for the sampling timepoint indicated.
* qpcr_dct_90k: deltaCt value for 90k qPCR, performed on RNA extracted from PBMCs. "no" indicates that this was not measured for the sampling timepoint indicated.
* who: The WHO classification of disease severity ascribed to a patient.
* sex: The sex of the patient. M = male, F = female.
* age: The age in years of a patient. 
* dex: Whether Dexamethasone was employed in the treatment of the patient, either "yes" or "no".

