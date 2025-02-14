#!/usr/bin/env python3
import argparse
import pysam
import sys
import gzip
import pprint
import time
from datetime import timedelta
import logging
from marktrimming import version
version = version.__version__
def main():
    marktrimming_cli()

def marktrimming_cli():
    start_time = time.time()
    # Configure logging
    logging.basicConfig(level=logging.INFO, format='# %(asctime)s ##INFO## %(message)s')
    logging.info("Started "+sys.argv[0]+" version "+version+".")
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        description='Takes queryname sorted reads and bam file and annotates the data with data from trimmed reads if present. ' + 
        'Should be used on an unaligned bam produced by fgbio or picard/gatk tools then the fastq dumped reads thoud be ran ' +
        'through the trimming tool or similar to produce the trimmed fq. This produces a bam similar to markIlluminaAdapters to be analysed with mergeBamAlignment. Tested on cutadapt without deleting reads.',
        epilog='Example useage python3 marktrimming.py --input input.bam --output output.bam --fastq qname_R1.fastq.gz --fastq qname_R2.fastq.gz')
    parser.add_argument('--output',help='Output bam. Should support "-" as /dev/stdout', required=True)
    parser.add_argument('--fastq', nargs='+' ,type=str, help='Input fastqs one or more arguments (two max). Should not be streamable.',action='append', required=True)
    parser.add_argument('--input', help='Input bam. Should support "-" as /dev/stdin', required=True )
    parser.add_argument('--ubam-out',help='Write uncompressed output bam. Helpful for speeding up on stream processing',  type=bool, required=False, default=False )
    parser.add_argument('--ubam-in',help='Read uncompressed input bam.', type=bool, required=False, default=False )
    parser.add_argument('--trimming-tag',help='Tag to add in trimmings. Default is XT for gatk compatibility', type=str, required=False, default="XT")
    parser.add_argument('--version', action='version', version='%(prog)s '+version)


    args = parser.parse_args()
    #main subroutine
    count,trimmed_count_first_of_pair,trimmed_count_second_of_pair = marktrimming(args)
    
    logging.info("Ended script.")
    logging.info(" Run time " + str(timedelta(seconds=time.time()- start_time)) + ".")
    logging.info(" Records processsed "+str(count) + ".")
    logging.info(" Trimmed first of pair "+ str(trimmed_count_first_of_pair) + ".")
    logging.info(" Trimmed second of pair "+ str(trimmed_count_second_of_pair) + ".")


def marktrimming(args):
    
    #pseudocode
    # handle switches
    # read trimmed fastq inputs keeping them (hopefully) synced by reading new records if their headers dont match with the bam data.
    # read initial fqs -> interate bam records -> sync fastq if header not eq to bam.
    # update XT tags if the read is shorter then the unaligned bam.  
    # assume bam format without @SQ reference sequence header.
    readmode = 'rb'
    inbam = pysam.AlignmentFile(args.input, readmode, check_sq=False, require_index=False)
    
    #exit(1)
    #there aint no logic to be found here
    header = pysam.samtools.view("-H",args.input)+ "\t".join(["@PG","ID:marktrimming", "PN:marktrimming", "VN:version"+version, "CL:"+' '.join(sys.argv)])
    
    writemode = 'wb'
    if args.ubam_out:
        writemode = 'wbu'
    outfile = pysam.AlignmentFile(args.output, writemode, text=header, check_sq=False, require_index=False)
    
    #intiatiate trimmed fastq reading and read them to an array of (paired) data. Every entry contains 4 lines of fastq data single line sequence/qual fastq data assumed. 
    fqs = []
    for list in args.fastq:
        for fqfile in list:
            #print("input fastqs "+ str(fqfile), file=sys.stderr)
            fqs.append( FastqFile(fqfile, 'r').open())

    count = 0
    fastqlineno = 0

    trimmed_count_first_of_pair = 0
    trimmed_count_second_of_pair = 0
    fastqrecords=[]
    for f in fqs:
        fastqrecords.append(f.read_record())
    fastqlineno = fastqlineno + 1
    #single bam record/end to multiple fastq records
    for lineno, bamrecord in enumerate(inbam):
        if lineno % 100000 == 0: 
            logging.info(" Records processsed "+str(lineno) + ".")
        def stripfqheader(header):
            #strips of the optional ["/1","/2", and space separated comments"
            header=header.lstrip('@')
            for char in [' ','/']:
                if char in header:
                    header=header[:header.index(char)]
            return header
        def checkpaired(records):
            #This checks of the fastq data is in sync structure like :  
            #[["@fq1 header","ATCGATCGN","+","BASEQUALS"],["@fq2 header","ATCGATCGN","+","BASEQUALS"]] 
            #checks if [*][0] is eq between array of records
            for record in records:
                if not stripfqheader(records[0][0]) == stripfqheader(record[0]):
                    print("error :: fastq files out of sync!! records:", file=sys.stderr)
                    print("R1"+ records[0][0] + "R2"+ record[1][0], file=sys.stderr)
                    print(records, file=sys.stderr)
                    exit(1)
            return True
        

        checkpaired(fastqrecords)
        #print(fastqrecords)
        
        def fastq_to_bam_iseq(fastq_records,bam_record):
            #This compares the header IDs to be the same between bam and fastq files
            #this might need a trycatch structure
            #try:
                if not bam_record.query_name == stripfqheader(fastq_records[0][0]):
                    #print("error :: fastq-bam file combi out of sync!! records:")
                    #print("error :: R1"+ stripfqheader(fastq_records[0][0]) + "error :: R2"+ bam_record.query_name)
                    return False
                return True
            #except:
            #    print("error :: fastq-bam file combi out of sync!! records:fq"+stripfqheader(fastq_records[0][0])+"bam"+bam_record.tostring())
            #    print("error :: ####Fastq\nerror :: ")
            #    print(fastq_records) 
            #    print ("\nerror :: ####Bam\nerror :: ")
            #    print( bam_record )
            #    exit(1) 

        #Sync fastq to bam bacause the fastq should be always lagging assuming the trimming might delete a read or so.
        sync_counter = 0
        sync_counter_max = 1000
        if not fastq_to_bam_iseq(fastqrecords,bamrecord):
            while (not fastq_to_bam_iseq(fastqrecords,bamrecord)) or (sync_counter > sync_counter_max):
                fastqrecords=[]
                for f in fqs:
                    fastqrecords.append(f.read_record())
                fastqlineno = fastqlineno + 1
                sync_counter = sync_counter + 1
                if sync_counter > sync_counter_max: 
                    print("error :: fastq-bam file combi out of sync!! records:", file=sys.stderr)
                    print("Fastq"+ (",".join([fastqrecords[0][0],fastqrecords[0][0]]))+"fastqlineno"+str(fastqlineno)  + "bam" + bamrecord.tostring()+"bamlineno" + str(lineno) , file=sys.stderr)
                    exit(1)
        #compare post trim with untrimmed lenght
        if  bamrecord.is_read1 and len(fastqrecords[0][1]) < len(bamrecord.query_sequence):
            #print(len(bamrecord.query_sequence), file=sys.stderr)
            #print(len(fastqrecords[0][1]), file=sys.stderr)
            #pprint.pprint(fastqrecords, stream=sys.stderr)
            bamrecord.set_tag(args.trimming_tag, len(fastqrecords[0][1])+1, value_type='i', replace=True)
            #print(bamrecord, file=sys.stderr)
            #exit(1)
            trimmed_count_first_of_pair = trimmed_count_first_of_pair + 1 
        if bamrecord.is_read2 and len(fastqrecords[1][1]) < len(bamrecord.query_sequence):
            #print(len(bamrecord.query_sequence), file=sys.stderr)
            #print(len(fastqrecords[1][1]), file=sys.stderr)
            #pprint.pprint(fastqrecords, stream=sys.stderr)
            bamrecord.set_tag(args.trimming_tag, len(fastqrecords[1][1])+1, value_type='i', replace=True)
            #print(bamrecord, file=sys.stderr)
            #exit(1)
            trimmed_count_second_of_pair = trimmed_count_second_of_pair + 1 

        outfile.write(bamrecord)
        count = count + 1
    #this is needed due to htslib needs to finish writing the last bam block.
    outfile.close()
    return count,trimmed_count_first_of_pair,trimmed_count_second_of_pair
    #notes
    test="""[umcg-mterpstra@gearshift Bestie]$ (ml Pysam ;python3 marktrimming.py --output tests/test.tmp.bam --input tests/runs/integration/cromwell-executions/FastqToVariants/4c1193d5-c447-4739-a5ff-bb1af87dc88e/call-fqToBam/shard-1/FastqToBam/*/call-fastqToUnmappedBam/shard-0/execution/*_ATTTTA+TAAAAT.L002_unaligned.bam --fastq tests/runs/integration/cromwell-executions/FastqToVariants/4c1193d5-c447-4739-a5ff-bb1af87dc88e/call-fqToBam/shard-1/FastqToBam/*/call-adaptertrim/shard-*/execution/MySample2_HXXXXXXXXXXXD_LL002_S02_trim_R1.fastq.gz --fastq tests/runs/integration/cromwell-executions/FastqToVariants/4c1193d5-c447-4739-a5ff-bb1af87dc88e/call-fqToBam/shard-1/FastqToBam/*/call-adaptertrim/shard-*/execution/MySample2_HXXXXXXXXXXXD_LL002_S02_trim_R2.fastq.gz )
[E::idx_find_and_load] Could not retrieve index file for 'tests/runs/integration/cromwell-executions/FastqToVariants/4c1193d5-c447-4739-a5ff-bb1af87dc88e/call-fqToBam/shard-1/FastqToBam/c54ccd72-af28-42f5-9387-c3aa26b804d8/call-fastqToUnmappedBam/shard-0/execution/002_HXXXXXXXXXXXD_ATTTTA+TAAAAT.L002_unaligned.bam'
tests/runs/integration/cromwell-executions/FastqToVariants/4c1193d5-c447-4739-a5ff-bb1af87dc88e/call-fqToBam/shard-0/FastqToBam/997a54a0-c1e3-44e0-aef1-681eceae6a73/call-adaptertrim/shard-0/execution/MySample_HXXXXXXXXXXXD_LL001_S01_trim_R1.fastq.gz.sorted.fastq.gz
r
tests/runs/integration/cromwell-executions/FastqToVariants/4c1193d5-c447-4739-a5ff-bb1af87dc88e/call-fqToBam/shard-0/FastqToBam/997a54a0-c1e3-44e0-aef1-681eceae6a73/call-adaptertrim/shard-0/execution/MySample_HXXXXXXXXXXXD_LL001_S01_trim_R2.fastq.gz.sorted.fastq.gz
r
[['@chr8:127736594-127740958_1000_1411_1:0:0_3:0:0_2b2cd/1', 'GCGCTCCCCCGGGAGCCCGGAGCGCAAAGCCCGGGAGTCGGCCGGGCAGCGGCAGAGG', '+', 'IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665'], ['@chr8:127736594-127740958_1000_1378_1:0:0_2:0:0_803d5/2', 'GGAGTCGGCCCCGCAGCGGCAGAGGAATCGCAATCGGCCCTGGCGCCGTTAAGAAGCC', '+', 'IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665']]
chr8:127736594-127740958_1000_1361_3:0:0_4:0:0_55db3    77      *       0       0       None    *       0       0       GAGGAATCGAAATCTGCCCTGGCGCCCTTAAGAAGGCGCGGGAGGTGGCGGTGAGGAAAACAATTTGGCAAAATCCAAGGCACAAAGTTTTGCGCCACCTGAAGGAGAAGGCGAGAGGCGCCTGG     array('B', [40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 37, 37, 37, 37, 37, 37, 37, 36, 36, 36, 36, 36, 36, 36, 35, 35, 35, 35, 35, 35, 34, 34, 34, 34, 34, 34, 33, 33, 33, 33, 33, 32, 32, 32, 32, 32, 31, 31, 31, 31, 30, 30, 30, 30, 30, 29, 29, 29, 29, 28, 28, 28, 28, 27, 27, 27, 27, 26, 26, 26, 25, 25, 25, 25, 24, 24, 24, 23, 23, 23, 23, 22, 22, 22, 21, 21, 21, 20, 20, 20])      [('RG', '002_HXXXXXXXXXXXD_ATTTTA+TAAAAT.L002')]

[umcg-mterpstra@gearshift Bestie]$ (ml SAMtools&& samtools view -h tests/runs/integration/cromwell-executions/FastqToVariants/4c1193d5-c447-4739-a5ff-bb1af87dc88e/call-fqToBam/shard-0/FastqToBam/997a54a0-c1e3-44e0-aef1-681eceae6a73/call-fastqToUnmappedBam/shard-0/execution/001_HXXXXXXXXXXXD_ATTTTA+TAAAAT.L001_unaligned.bam| head; ) 
@HD     VN:1.6  SO:queryname
@RG     ID:001_HXXXXXXXXXXXD_ATTTTA+TAAAAT.L001 SM:MySample     LB:MySampleATTTTA+TAAAAT        PL:ILLUMINA     PU:001_ATTTTA+TAAAAT.L001       DT:2024-03-12T00:00:00+0000
@PG     ID:samtools     PN:samtools     VN:1.19.2       CL:samtools view -h tests/runs/integration/cromwell-executions/FastqToVariants/4c1193d5-c447-4739-a5ff-bb1af87dc88e/call-fqToBam/shard-0/FastqToBam/997a54a0-c1e3-44e0-aef1-681eceae6a73/call-fastqToUnmappedBam/shard-0/execution/001_HXXXXXXXXXXXD_ATTTTA+TAAAAT.L001_unaligned.bam
chr8:127736594-127740958_1000_1378_1:0:0_2:0:0_803d5    77      *       0       0       *       *       0       0       GTGTTAAAGCCCGCGGCTGAGCTCGCCACTCCAGCCGGCGAGAGAAAGAAGCAAAGCT      IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665        RG:Z:001_HXXXXXXXXXXXD_ATTTTA+TAAAAT.L001
chr8:127736594-127740958_1000_1378_1:0:0_2:0:0_803d5    141     *       0       0       *       *       0       0       GGAGTCGGCCCCGCAGCGGCAGAGGAATCGCAATCGGCCCTGGCGCCGTTAAGAAGCC      IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665        RG:Z:001_HXXXXXXXXXXXD_ATTTTA+TAAAAT.L001
chr8:127736594-127740958_1000_1380_2:0:0_0:0:0_999b6    77      *       0       0       *       *       0       0       GTGTTAAAGCCCGCGTCTGAGCTCGCCAGTCCAGCCGGCGAGAGAAAGAAGAAAAGCT      IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665        RG:Z:001_HXXXXXXXXXXXD_ATTTTA+TAAAAT.L001
chr8:127736594-127740958_1000_1380_2:0:0_0:0:0_999b6    141     *       0       0       *       *       0       0       CGGGAGTCGGCCCCGCAGCGGCAGAGGAATCGAAATCGGCCCTGGCGCCCTTAAGA        IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::98876  RG:Z:001_HXXXXXXXXXXXD_ATTTTA+TAAAAT.L001
chr8:127736594-127740958_1000_1389_2:0:0_1:0:0_7e197    77      *       0       0       *       *       0       0       CGCAAAGCCCGGGAGTCGGCCCCGCAGCGGCAGAGGAATCGAAATCTGCCCTGGCGCC      IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665        RG:Z:001_HXXXXXXXXXXXD_ATTTTA+TAAAAT.L001
chr8:127736594-127740958_1000_1389_2:0:0_1:0:0_7e197    141     *       0       0       *       *       0       0       GTGTTAAAGCCCGCGGGTGAGCTCGCCACTCCAGCCGGCGAGAGAAAGAAGAAACGCT      IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665        RG:Z:001_HXXXXXXXXXXXD_ATTTTA+TAAAAT.L001
chr8:127736594-127740958_1000_1391_1:0:0_2:0:0_e5ae1    77      *       0       0       *       *       0       0       GTGTTAAAGCCCGCGGCTGAGCTCGCCACTCCAGCCGGCGAGAGAAAGAAGAAAAGC       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::988766 RG:Z:001_HXXXXXXXXXXXD_ATTTTA+TAAAAT.L001

[umcg-mterpstra@gearshift Bestie]$ (ml HTSlib&& bgzip -dc tests/runs/integration/cromwell-executions/FastqToVariants/4c1193d5-c447-4739-a5ff-bb1af87dc88e/call-fqToBam/shard-0/FastqToBam/
997a54a0-c1e3-44e0-aef1-681eceae6a73/call-adaptertrim/shard-0/execution/MySample_HXXXXXXXXXXXD_LL001_S01_trim_R1.fastq.gz.sorted.fastq.gz | paste - - - - | head
> ) 
@chr8:127736594-127740958_1000_1378_1:0:0_2:0:0_803d5/1 GTGTTAAAGCCCGCGGCTGAGCTCGCCACTCCAGCCGGCGAGAGAAAGAAGCAAAGCT      +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665
@chr8:127736594-127740958_1000_1380_2:0:0_0:0:0_999b6/1 GTGTTAAAGCCCGCGTCTGAGCTCGCCAGTCCAGCCGGCGAGAGAAAGAAGAAAAGCT      +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665
@chr8:127736594-127740958_1000_1389_2:0:0_1:0:0_7e197/1 CGCAAAGCCCGGGAGTCGGCCCCGCAGCGGCAGAGGAATCGAAATCTGCCCTGGCGCC      +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665
@chr8:127736594-127740958_1000_1391_1:0:0_2:0:0_e5ae1/1 GTGTTAAAGCCCGCGGCTGAGCTCGCCACTCCAGCCGGCGAGAGAAAGAAGAAAAGC       +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::988766
@chr8:127736594-127740958_1000_1395_2:0:0_0:0:0_b342b/1 CCGGAGCGCAAAGCCCGGGAGTCGGCCCCGCAGCGGCAGAGGAATCGAAATCGGCCCT      +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665
@chr8:127736594-127740958_1000_1401_0:0:0_0:0:0_e339c/1 GGGAGCCCGGAGCGCAAAGCCCGGGAGTCGGCCCCGCAGCGGCAGAGGAATCGAAATC      +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665
@chr8:127736594-127740958_1000_1401_0:0:0_1:0:0_800db/1 GGGAGCCCGGAGCGCAAAGCCCTGGAGTCGGCCCCGCAGCGGCAGAGGAATCGAAATC      +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665
@chr8:127736594-127740958_1000_1403_1:0:0_0:0:0_58712/1 CCGGGAGCCCGGAGCGCAAAGCCCGGGAGTCGGCCCCGCAGCGGCAGAGGAATCGAA       +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::988766
@chr8:127736594-127740958_1000_1405_1:0:0_2:0:0_be547/1 GTGTTAAAGCCCGCGGCTGAGCTCGCCAGTCCAGCCGGCGAGAGAAAGAAGAAAAGCT      +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665
@chr8:127736594-127740958_1000_1406_1:0:0_1:0:0_ef320/1 GTGTTAAAGCCCGCGGCTGCGCTCGCCACTCCAGCCGGCGAGAGAAAGAAGAAAAGCT      +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665
[umcg-mterpstra@gearshift Bestie]$ (ml HTSlib&& bgzip -dc tests/runs/integration/cromwell-executions/FastqToVariants/4c1193d5-c447-4739-a5ff-bb1af87dc88e/call-fqToBam/shard-0/FastqToBam/997a54a0-c1e3-44e0-aef1-681eceae6a73/call-adaptertrim/shard-0/execution/MySample_HXXXXXXXXXXXD_LL001_S01_trim_R2.fastq.gz.sorted.fastq.gz | paste - - - - | head; ) 
@chr8:127736594-127740958_1000_1378_1:0:0_2:0:0_803d5/2 GGAGTCGGCCCCGCAGCGGCAGAGGAATCGCAATCGGCCCTGGCGCCGTTAAGAAGCC      +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665
@chr8:127736594-127740958_1000_1380_2:0:0_0:0:0_999b6/2 CGGGAGTCGGCCCCGCAGCGGCAGAGGAATCGAAATCGGCCCTGGCGCCCTTAAGA        +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::98876
@chr8:127736594-127740958_1000_1389_2:0:0_1:0:0_7e197/2 GTGTTAAAGCCCGCGGGTGAGCTCGCCACTCCAGCCGGCGAGAGAAAGAAGAAACGCT      +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665
@chr8:127736594-127740958_1000_1391_1:0:0_2:0:0_e5ae1/2 AGCGCAAAGCCCGGGAGTCGGCCCCGCAGCGGCAGAGGACTCGAAATCTGCCCTGGCG      +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665
@chr8:127736594-127740958_1000_1395_2:0:0_0:0:0_b342b/2 GTGTTAAAGCCCGCGGCTGAGCTCGCCACACCAGCCGGCGAGAGAAAGAAGAAAAGCT      +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665
@chr8:127736594-127740958_1000_1401_0:0:0_0:0:0_e339c/2 GTGTTAAAGCCCGCGGCTGAGCTCGCCACTCCAGCCGGCGAGAGAAAGAAGAAAAGCT      +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665
@chr8:127736594-127740958_1000_1401_0:0:0_1:0:0_800db/2 GTGTTAAAGCCCGCGGCTGAGCTCGCCACTCCAGCCGGCGAGAGAAAGAAGAAAAGCT      +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665
@chr8:127736594-127740958_1000_1403_1:0:0_0:0:0_58712/2 GTGTTAAAGCCCGCGGCTGAGCTCGCCACTCCAGCCGGCGAGAGAAAGAAGAAAAGCT      +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665
@chr8:127736594-127740958_1000_1405_1:0:0_2:0:0_be547/2 CCCCGGGAGCCCGGCGCGCAAAGCCCGGGAGTCGGCCCCGCAGCGGCAGCGGAATCG       +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::988766
@chr8:127736594-127740958_1000_1406_1:0:0_1:0:0_ef320/2 CCCCCGGGAGCCCGGAGCGCAAAGCCCGGGAGTCGGCCCCGCCGCGGCAGAGGAATCG      +       IIIIIIIIIIHHHHHHGGGGGFFFEEEEDDCCCBBBAA@@??>>==<;;::9887665

"""

class FastqFile:
    def __init__(self, filename, mode='r'):
        self.filename = filename
        self.mode = mode
        self.file_handle = None
        self.buffered_record = None
    
    def open(self):
        #print(self.filename)
        #print(self.mode)
        self.file_handle = None
        if self.filename.endswith('.gz'):
            self.file_handle = gzip.open(self.filename, self.mode + 't')
        else:
            self.file_handle = open(self.filename, self.mode)
        
        return(self)
    
    def close(self):
        if self.file_handle is not None:
            self.file_handle.close()
            self.file_handle = None
    
    def read_record(self):
        if self.buffered_record:
            record = self.buffered_record
            self.buffered_record = None
            return record
        
        if self.mode != 'r':
            raise ValueError("File not opened in read mode")
        lines = [self.file_handle.readline() for _ in range(4)]
        for line in lines[1:]:
            if len(line) == 0:
                raise ValueError("incorrect or corrupt fastq format.\nrecord:"+"\n".join(lines))
        if len(lines[0]) == 0:
            raise ValueError("Trying to read beyond eof: time for some more sanity checks!\nrecord:"+"\n".join(lines))
            return ["","","",""]
        lines = [_.rstrip() for _ in lines]
        return lines
    
    def unread_record(self, record):
        if self.buffered_record:
            raise ValueError("Cannot buffer more than one record")
        self.buffered_record = record
    
    def write_record(self, record):
        if self.mode != 'w':
            raise ValueError("File not opened in write mode")
        if len(record) != 4:
            raise ValueError("Invalid record format")
        for line in record:
            self.file_handle.write(line + '\n')
    
    def __enter__(self):
        self.open()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


# __name__
if __name__=="__main__":
    main()
    
