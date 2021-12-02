import sys, os

in_file = sys.argv[1]
with open(in_file, 'r') as file:
    for line in file:
        line = line.strip()
        if '##fileformat' in line:
            print('##fileformat=VCFv4.3')
        elif '##INFO=<ID=AF' in line:
            print(
                '##INFO=<ID=AF,Number=A,Type=Float,'
                'Description="Allele Frequency">')
        elif '#CHROM' in line:
            print(
                '##contig=<ID=NC_045512.2,length=29903,'
                'URL=https://www.ebi.ac.uk/ena/browser/api/fasta/MN908947.3>')
            print(line)
        else:
            print(line)
