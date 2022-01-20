# restriction_site_generator
Site Specific Restriction Site Generator

********************************************************************



*************************************************************************
Restriction site Generator is a tool which allows the creation of restriction site along with or without mutation at specific site in primers. 
Applications: Creation of restriction site in primers to detect the mutation, CRISPR Cas insertions or replacements for downstream analysis with restriction enzyme digestion. This is alternative to PCR detection method. 

*********************************************************************
COMMANDLINE interface:


Requirements: 
Biopython,
Numpy,
Pandas


STEPS:
Step1: Enter the primer sequence (max 100nt)
Step2: Choose the correct in-frame ORF
Step3: Choose Amino acid to mutate by entering the amino acid number in displayed sequence 
Step4: Choose the amino acid only from +3 to -3; the algorithm requires two amino acids before and after the selected amino acid to maximize the restriction sites. If not enough sequence, please increase primer length. E.g., “GNGHHGHGHGHGSG,” to create mutation and restriction site at “N” or “S” will need one more amino acid at N-terminal or C-terminal.
Results will be created in excel file “results.xlsx” downloaded same folder. 
