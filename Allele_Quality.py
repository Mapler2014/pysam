#!/usr/bin/env python
import pysam
import sys,argparse

usage='''
 Python Allele_Quality.py --bam=sample.bam  --tsv=chr.tsv  --out=allele_summary.txt
'''

def paramsParse():
    parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--bam', help = "bam file", required = True)
    parser.add_argument('--tsv', help = "Two column tab seperated file that contains the query locations", required = True)
    parser.add_argument('--out', help = "output file", required = True)
    return parser.parse_args()

#Pysam was used to read and manipulate mapped short read sequence data stored in SAM/BAM files
def allele(bam,chrID,position):
    samfile=pysam.AlignmentFile(bam,"rb")
    chrid=chrID
    chr_pos=int(position)
    start_pos=chr_pos-50
    end_pos=chr_pos+50
    allele_dict=dict()
    allele_summary=list()
    for pileupcolumn in samfile.pileup(chrid,start_pos,end_pos,min_base_quality=0):
        if(pileupcolumn.pos==chr_pos-1):
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    # query position is None if is_del or is_refskip is set.
                    q=pileupread.alignment.query_qualities
                    allele=pileupread.alignment.query_sequence[pileupread.query_position]
                    quality=q[pileupread.query_position]
                    if allele in allele_dict:
                        allele_dict[allele].append(quality)
                    else:
                        allele_dict[allele]=[quality]
                   
                else:
                    print(pileupread.alignment.query_name,"is indel\n")
    samfile.close()
    for k in allele_dict.keys():
        qua_list=allele_dict[k]
        ave_qua=format(sum(qua_list)/len(qua_list),".2f")
        allele_summary.extend([k,str(ave_qua)])
    return allele_summary


if __name__=='__main__':
    args=paramsParse()
    FH=open(args.tsv,'r')
    OUT=open(args.out,'w')
    # The header in the output need be modified by the columns of output if there are more than two alleles
    OUT.write("Chromosome\tBase\tAllele1\tAvg_1\tAllele2\tAvg_2\n")
    for line in FH:
        line_list=line.strip().split("\t")
        if line_list[0] != "Chromosome":
            allele_sum=allele(args.bam, line_list[0],line_list[1])
            out_string=line.strip()+"\t"+"\t".join(allele_sum)+"\n"
            OUT.write(out_string)
    FH.close()
    OUT.close()



