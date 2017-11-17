***********
DESCRIPTION
***********
Smith Waterman algorithm for Local Sequence Alignment. Computes alignment with both affine gap penalty and naive gap penalty.
Affine Gap Penalty: The gap opening penalty is 16. Gap Extension penalty is 4.
Navice Gap Penalty: Gap penalty is 4.

*****
Usage
*****
1) The directory must contain an 'input.fa' file in the fasta file format. This name is hardcoded in the code.
2) Once the program has been compiled, execute the bash script 'output.sh'. It parses the output and puts it in appropriate files.

********************
Notes On Compilation
********************
The program uses linux system call, therefore it might not work in a windows environment.
