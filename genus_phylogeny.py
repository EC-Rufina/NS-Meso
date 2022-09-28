import Bio
from Bio import SeqIO
from Bio import SeqFeature
import sys
import os
import csv
import random

## Genus phylogeny
#this script needs to be modified

#before starting you need to create 2 folders:
#one folder called proteinortho, where you add the proteinortho result .tsv file
#one folder called prokka_annotation where you place all genomes annotation from prokka (we just need .faa and .gbk for this)


##before you start type below the name of your .tsv file#
name_of_pr_ortho_file="proteinortho_file.tsv"
input_path = "proteinortho/"

#then change the following paramenter: how many genomes are you using?
z = 190

#name of the path to the prokka annotation
input_path_prokka="prokka_input/"


#now you don't need to change anything else#
#create all the folders
output_path_RAW= "coregenes_locustags/"
final_ouput_path = "coregenes/"

os.mkdir(output_path_RAW)
os.mkdir(final_ouput_path)
os.mkdir("coregenes_alignment_prank")
os.mkdir("coregenes_alignment_prank_cleaned")

proteinortho_results = csv.reader(open(input_path + name_of_pr_ortho_file, 'r'), delimiter='\t')
x = 3
y = x + z

r = range(3, y)


#now we need to create a list with the names of the genomes used in proteinortho
ICElist=[]
for row in proteinortho_results:
    for i in range(x, y):
        
        ICElist.append(row[i])
    break

ICElist = [s.replace(".faa", "") for s in ICElist]


#we extract the locus tag of each single-copy conserved gene
for line in proteinortho_results:
    
    if int(line[0])==z and int(line[1])==z:
        #print('extracting locus tags from	 ' + line[3])
        out1 = open(output_path_RAW + line[3] + '_singlegene.txt', 'w')
        for i in r:
            out1.write(line[x] + ',')
            x +=1
        x = 3 
        out1.close()    


#yes we did it already
proteinortho_results = csv.reader(open(input_path + name_of_pr_ortho_file, 'r'), delimiter='\t')
header = next(proteinortho_results)
	
print('all locus tags of the core genes have been identified, now I will extract the sequences' )

num = 1


for nome in os.listdir(output_path_RAW):
    ofile = open(final_ouput_path + nome + "_coregenes.fasta", "w")

    mylist=open(output_path_RAW + nome , 'r')
    coregenes = mylist.read().split(',')
    
    print('I am now extracting from all the genomes the gene number   ' + str(num))
    num +=1
    for element in ICElist:
        #print('parsing   ' + element)
               
        ofile.write(">" + str(element) + "\n")
    

###this loop here is made so the genes are in the same order, does not follow the order the genes are found on the genome###
        
        gbk_input = SeqIO.parse(input_path_prokka + element + ".gbk", "genbank")
        #print("i am parsing  " + element)
        for genome in gbk_input:
            for gene in genome.features:
                if(gene.type =="CDS"):
                    locustag=gene.qualifiers['locus_tag'][0]
                    prseq=gene.qualifiers['translation'][0]
                    DNAseq=gene.extract(genome.seq)
                    if locustag in coregenes:
                        ofile.write(str(DNAseq) + "\n")
                        continue
                    continue
                continue
            continue
        continue

    ofile.close()




#now we align with prank all the core genes

for filename in os.listdir(final_ouput_path):
   
    if filename.endswith(".fasta"):
        name=filename.replace('.fasta', '')
        #command='mafft --thread 16 --auto {}{} > coregenes_alignment/al_{}.fasta'.format(final_ouput_path,filename,name)
        #os.system(command)
        command_prank='prank -d={}{} -o=coregenes_alignment_prank/al_{}.fasta'.format(final_ouput_path,filename,name)
        os.system(command_prank)

#we then remove gaps


working_directory='coregenes_alignment_prank/'
for filename in os.listdir(working_directory):
    name=filename.replace('.fasta', '')
    print("cleaning    " + name)
    command='goalign clean sites -i {}{} > coregenes_alignment_prank_cleaned/{}.fasta'.format(working_directory,filename,name)
    os.system(command)




#we concatenate the alignments 


command='goalign concat -i coregenes_alignment_prank_cleaned/*.fasta -o PRANK_clean_concatenated_core_aln.fasta'
os.system(command)
#we finally make the tree

random1=random.randrange(0,32767)
random2=random.randrange(0,32767)

###here you want to change the -o tag, that is your ougroup strain 
#command='raxmlHPC -f a -p {} -x {} -N 100 -m GTRCATX -T 16 -s clean_concatenated_core_aln.fasta -n raxml_comp_cleaned_coregenes -o Rhizobium_leguminosarum_bv_viciae_USDA2370'.format(random1,random2) 
#os.system(command)
command='raxmlHPC -f a -p {} -x {} -N 100 -m GTRCATX -T 16 -s PRANK_clean_concatenated_core_aln.fasta -n PRANK_raxml_comp_cleaned_coregenes -o M.sp.--HJP'.format(random1,random2) 
os.system(command)


