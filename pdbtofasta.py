import sys

"""
The purpose of the script is to extract the residue sequence from the atomic coordinates of a protein.
If the protein has multiple chains, the fasta output separates the chains by "/". 
The description line has ">" in the first coloumn followed by th total number of residues 
in the PDB file. The output is set to contain a maximum of 70 characters in each line.
"""

if len(sys.argv) <= 1 :
    print 'USAGE : python pdbtofasta.py input.pdb > output.fasta'
    exit ()
    
input_file = open(sys.argv[1])
input_pdb = input_file.readlines()

res_dict = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
           'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
           'TYR':'Y','VAL':'V'}

prev_res = -1
fasta = ""
total_res = 0
for line in input_pdb:
    if (line[0:4] == "ATOM" and len(line)>53):
        res_number = eval(line[22:26])
        res_type = line[17:20]

        if (res_number > prev_res):
            fasta = fasta+res_dict[res_type]
            total_res+=1
            if total_res%70==0:
                fasta=(fasta+"\n")
        prev_res = res_number

    if (line[0:3]=="TER"):
        fasta = fasta + '/' 
        prev_res = -1 

print ('>%s %d' %(sys.argv[1],total_res))
print (fasta)



