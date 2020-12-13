#This code is for converting gene into RNA,
#             for splicing exons out of RNA,
#             for converting exons into protein.

# convert gene into RNA
gene = 'ACATACTCGTTATTTACTTGGTTGCTACAGGTCCGGGACAGGCAAGATCCTTTGAACACCGTCGAGTAGTGACGCGTAGTGGGGGACCTTAAAGGGCAGACAAGCTATTCTGCAAAGAGAATAACAGGAAGGTCTCACGATAGTCTTTGTTCAATCCCGTTTGAAGCTCTCTCACTGCTCCAGGTTAAGACCGGTTGTTCCCAACTACGAACGAGGAGTTGTCCAAATCGTACTGATGAGCCTTTTCACTCTTGCAGGATAGAGCTACGCTTCATTGAGTGAGACAGAGATTCGAGCTCAGCCGTCGTGCTCACCGCTACAATCAAACGTTATTCCGTGGACCCAGAATGGGCCGTCAATGGGAACCAGTCGTGGGAATGCGGGGGGCGGAGGACCGCTA'
# Enter your own sequence text file destination between ' and '
replace_a=(gene.replace('A','u'))
replace_t=(replace_a.replace('T','a'))
replace_c=(replace_t.replace('C','g'))
replace_g=(replace_c.replace('G','c'))
rna = (replace_g.upper())
print('RNA CODE','=',rna)

# splice exons out of RNA
exon1 = rna[5:58]
exon2 = rna[72:133]
exon3 = rna[190:276]
exon4 = rna[340:398]
exons = exon1 + exon2 + exon3 + exon4
print('SPLICED RNA','=',exons)

# convert exons into protein
gencode = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
def translate_gene(exons):
    last_codon_start = len(exons) - 2
    protein = ""
    for start in range(0,last_codon_start,3):
        codon = exons[start:start+3]
        aa = gencode.get(codon.upper(), 'X')
        protein = protein + aa
    return protein
    
print('PROTEIN CODE','=',translate_gene(exons))