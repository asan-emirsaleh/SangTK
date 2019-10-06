# SangTK

## Dependencies
* biopython
* pandas
* sklearn

## Usage
```
python sang.py [-h]
```
#### Essential options
**```-d/--ab1_directory```**   
Directory containing .ab1 files to be converted into fasta file.   
   
**```-f/--ab1_file```**   
Single .ab1 file to be converted into fasta file.   
   
#### Optional options
**```-s/--split```**   
Split output into separate fasta files. Input as true or false. 
   
**```-o/--fa_name```**   
Input non-default name for fasta file (not valid if inputting >1 .ab1 file with -s flag).  
   
**```-pn/predict_nucleotide```**   
Use predictive algorithm to determine sequence. Calls nucleotides given peaks in .ab1 file. 

**```-p/predict_peak_and_nucleotide```**   
Use predictive algorithm to determine sequence. Calls both peaks and nucleotides. 
