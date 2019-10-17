# SangTK

## Function
SangTK is a toolkit for AB1 (sanger sequencing) files that can perform any combination of the following functions:
* Convert .abi files to .fa files
* Convert a directory of .abi files into a single fa file or separate .fa files
* Perform logistic regression on nucleotide peaks to improve sequence quality
* Call peaks using a Bi-Directional LSTM RNN
* Improve the quality of the supplied peak calls

## Dependencies
* python 3
* biopython
* pandas
* sklearn
* tensorflow/keras
* scipy
* numpy

## Usage
```
python sang.py [-h]
```
#### Essential (exclusive) inputs
**```-d/--ab1_directory```**   
Directory containing .ab1 files to be converted into fasta file.   
   
**```-f/--ab1_file```**   
Single .ab1 file to be converted into fasta file.   
   
#### Optional inputs
**```-s/--split```**   
Split output into separate fasta files. Input as true or false. 
   
**```-o/--fa_name```**   
Input non-default name for fasta file (not valid if inputting >1 .ab1 file with -s flag).  
   
**```-pn/--predict_nucleotide```**   
Use predictive algorithm to determine sequence. Calls nucleotides given peaks in .ab1 file. 

**```-p/--predict_peak_and_nucleotide```**   
Use predictive algorithm to determine sequence. Calls both peaks and nucleotides. 
