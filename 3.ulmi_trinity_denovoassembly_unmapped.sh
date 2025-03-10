#!/bin/bash
#-------------------------------------------------------------------------------
##############################################################
#                                                            #
#            Script for preprocessing rawdata                #
#                                                            #
##############################################################
#-------------------------------------------------------------------------------
# https://github.com/vchano/
# Licence: GNU General Public Licence Version 3
#-------------------------------------------------------------------------------

# STRUCTURE DATA BY GENOTYPE

# MDV1 EXP
mdv1e='Gen-004 Gen-005 Gen-007 Gen-012 Gen-016 Gen-021 Gen-027 Gen-034 Gen-036 Gen-038 Gen-049 Gen-061 Gen-063 Gen-065 Gen-066 Gen-081 Gen-103 Gen-104 Gen-107 Gen-113 Gen-117 Gen-119 Gen-134 Gen-135 Gen-139 Gen-152 Gen-157 Gen-163 Gen-170 Gen-179 Gen-180 Gen-193'
for sample in ${mdv1e}
do
mkdir $GENESIS/ALIGN/MDV1/${sample}
cp -a $GENESIS/MAPPING/EXP/${sample}/contam/fungi/un-conc-mate.1.fastq.gz $GENESIS/ALIGN/MDV1/${sample}
cp -a $GENESIS/MAPPING/EXP/${sample}/contam/fungi/un-conc-mate.2.fastq.gz $GENESIS/ALIGN/MDV1/${sample}
done

# MDV1 CONTROL
mdv1c='Gen-010 Gen-011 Gen-019 Gen-025 Gen-042 Gen-045 Gen-047 Gen-051 Gen-052 Gen-054 Gen-062 Gen-076 Gen-078 Gen-079 Gen-087 Gen-088 Gen-089 Gen-093 Gen-109 Gen-110 Gen-121 Gen-132 Gen-147 Gen-151 Gen-156 Gen-158 Gen-166 Gen-167 Gen-168 Gen-169 Gen-177 Gen-183'
for sample2 in ${mdv1c}
do
mkdir $GENESIS/ALIGN/MDV1/${sample2}
cp -a $GENESIS/MAPPING/CONTROL/${sample2}/fungi/un-conc-mate.1.fastq.gz $GENESIS/ALIGN/MDV1/${sample2}
cp -a $GENESIS/MAPPING/CONTROL/${sample2}/fungi/un-conc-mate.2.fastq.gz $GENESIS/ALIGN/MDV1/${sample2}
done


# MDV23 EXP
mdv23e='Gen-008 Gen-022 Gen-026 Gen-028 Gen-031 Gen-035 Gen-037 Gen-044 Gen-048 Gen-067 Gen-071 Gen-073 Gen-075 Gen-080 Gen-085 Gen-086 Gen-090 Gen-091 Gen-095 Gen-118 Gen-123 Gen-142 Gen-148 Gen-154 Gen-155 Gen-159 Gen-160 Gen-165 Gen-178 Gen-184 Gen-185 Gen-191'
for sample3 in ${mdv23e}
do
mkdir $GENESIS/ALIGN/MDV23/${sample3}
cp -a $GENESIS/MAPPING/EXP/${sample3}/contam/fungi/un-conc-mate.1.fastq.gz $GENESIS/ALIGN/MDV23/${sample3}
cp -a $GENESIS/MAPPING/EXP/${sample3}/contam/fungi/un-conc-mate.2.fastq.gz $GENESIS/ALIGN/MDV23/${sample3}
done

# MDV23 CONTROL
mdv23c='Gen-009 Gen-033 Gen-046 Gen-050 Gen-056 Gen-057 Gen-064 Gen-068 Gen-069 Gen-074 Gen-082 Gen-084 Gen-092 Gen-096 Gen-100 Gen-102 Gen-105 Gen-106 Gen-112 Gen-115 Gen-130 Gen-133 Gen-149 Gen-153 Gen-164 Gen-171 Gen-176 Gen-181 Gen-186 Gen-188 Gen-189 Gen-190'
for sample4 in ${mdv23c}
do
mkdir $GENESIS/ALIGN/MDV23/${sample4}
cp -a $GENESIS/MAPPING/CONTROL/${sample4}/fungi/un-conc-mate.1.fastq.gz $GENESIS/ALIGN/MDV23/${sample4}
cp -a $GENESIS/MAPPING/CONTROL/${sample4}/fungi/un-conc-mate.2.fastq.gz $GENESIS/ALIGN/MDV23/${sample4}
done


# VAD2 EXP
vad2e='Gen-001 Gen-003 Gen-006 Gen-014 Gen-024 Gen-029 Gen-032 Gen-055 Gen-058 Gen-059 Gen-060 Gen-094 Gen-098 Gen-101 Gen-114 Gen-116 Gen-120 Gen-122 Gen-124 Gen-125 Gen-126 Gen-127 Gen-128 Gen-136 Gen-137 Gen-140 Gen-141 Gen-145 Gen-150 Gen-175 Gen-182 Gen-194'
for sample5 in ${vad2e}
do
mkdir $GENESIS/ALIGN/VAD2/${sample5}
cp -a $GENESIS/MAPPING/EXP/${sample5}/contam/fungi/un-conc-mate.1.fastq.gz $GENESIS/ALIGN/VAD2/${sample5}
cp -a $GENESIS/MAPPING/EXP/${sample5}/contam/fungi/un-conc-mate.2.fastq.gz $GENESIS/ALIGN/VAD2/${sample5}
done

# VAD2 CONTROL
vad2c='Gen-002 Gen-013 Gen-015 Gen-017 Gen-018 Gen-020 Gen-023 Gen-030 Gen-039 Gen-040 Gen-041 Gen-043 Gen-070 Gen-077 Gen-083 Gen-097 Gen-099 Gen-108 Gen-111 Gen-129 Gen-131 Gen-138 Gen-143 Gen-144 Gen-146 Gen-161 Gen-172 Gen-173 Gen-174 Gen-187 Gen-192 Gen-195'
for sample6 in ${vad2c}
do
mkdir $GENESIS/ALIGN/VAD2/${sample6}
cp -a $GENESIS/MAPPING/CONTROL/${sample6}/fungi/un-conc-mate.1.fastq.gz $GENESIS/ALIGN/VAD2/${sample6}
cp -a $GENESIS/MAPPING/CONTROL/${sample6}/fungi/un-conc-mate.2.fastq.gz $GENESIS/ALIGN/VAD2/${sample6}
done

# PREPARE DATA FILES AND RUN KMERGENIE
gunzip -r $GENESIS/ALIGN/MDV1

find $GENESIS/ALIGN/MDV1 -type f -name un-conc-mate.1.fastq -exec cat {} + > $GENESIS/ALIGN/MDV1/mdv1_1.fastq
gzip $GENESIS/ALIGN/MDV1/mdv1_1.fastq
find $GENESIS/ALIGN/MDV1/ -name un-conc-mate.1.fastq -type f -delete

kmergenie $GENESIS/ALIGN/MDV1/mdv1_1.fastq.gz -o mdv1_1 -t 8

find $GENESIS/ALIGN/MDV1 -type f -name un-conc-mate.2.fastq -exec cat {} + > $GENESIS/ALIGN/MDV1/mdv1_2.fastq
gzip $GENESIS/ALIGN/MDV1/mdv1_2.fastq
find $GENESIS/ALIGN/MDV1/ -name un-conc-mate.2.fastq -type f -delete
kmergenie $GENESIS/ALIGN/MDV1/mdv1_2.fastq.gz -o mdv1_2 -t 8


gunzip -r $GENESIS/ALIGN/MDV23

find $GENESIS/ALIGN/MDV23 -type f -name un-conc-mate.1.fastq -exec cat {} + > $GENESIS/ALIGN/MDV23/mdv23_1.fastq
gzip $GENESIS/ALIGN/MDV23/mdv23_1.fastq
find $GENESIS/ALIGN/MDV23/ -name un-conc-mate.1.fastq -type f -delete

kmergenie $GENESIS/ALIGN/MDV23/mdv23_1.fastq.gz -o mdv23_1 -t 8

find $GENESIS/ALIGN/MDV23 -type f -name un-conc-mate.2.fastq -exec cat {} + > $GENESIS/ALIGN/MDV23/mdv23_2.fastq
gzip $GENESIS/ALIGN/MDV23/mdv23_2.fastq
find $GENESIS/ALIGN/MDV23/ -name un-conc-mate.2.fastq -type f -delete
kmergenie $GENESIS/ALIGN/MDV23/mdv23_2.fastq.gz -o mdv23_2 -t 8


gunzip -r $GENESIS/ALIGN/VAD2

find $GENESIS/ALIGN/VAD2 -type f -name un-conc-mate.1.fastq -exec cat {} + > $GENESIS/ALIGN/VAD2/vad2_1.fastq
gzip $GENESIS/ALIGN/VAD2/vad2_1.fastq
find $GENESIS/ALIGN/VAD2/ -name un-conc-mate.1.fastq -type f -delete

kmergenie $GENESIS/ALIGN/VAD2/vad2_1.fastq.gz -o vad2_1 -t 8

find $GENESIS/ALIGN/VAD2 -type f -name un-conc-mate.2.fastq -exec cat {} + > $GENESIS/ALIGN/VAD2/vad2_2.fastq
gzip $GENESIS/ALIGN/VAD2/vad2_2.fastq
find $GENESIS/ALIGN/VAD2/ -name un-conc-mate.2.fastq -type f -delete
kmergenie $GENESIS/ALIGN/VAD2/vad2_2.fastq.gz -o vad2_2 -t 8

# RUNG TRINITY BY GENOTYPE

Trinity --seqType fq --max_memory 64G --CPU 8 --left $GENESIS/ALIGN/MDV1/mdv1_1.fastq.gz --right $GENESIS/ALIGN/MDV1/mdv1_2.fastq.gz --KMER_SIZE 32 --output $GENESIS/ALIGN/MDV1/TRINITY32/ 
Trinity --seqType fq --max_memory 64G --CPU 8 --left $GENESIS/ALIGN/MDV23/mdv23_1.fastq.gz --right $GENESIS/ALIGN/MDV23/mdv23_2.fastq.gz --KMER_SIZE 32 --output $GENESIS/ALIGN/MDV23/TRINITY32/ 
Trinity --seqType fq --max_memory 64G --CPU 8 --left $GENESIS/ALIGN/VAD2/vad2_1.fastq.gz --right $GENESIS/ALIGN/VAD2/vad2_2.fastq.gz --KMER_SIZE 32 --output $GENESIS/ALIGN/VAD2/TRINITY32/ 

# DOWNLOAD 454-PYROSEQUENCING (PERDIGUERO ET AL., 2015)

list='SRR1687227 SRR514281 SRR514280 SRR514277 SRR514276 SRR514275 SRR514274 SRR514253 SRR513995'
for item in ${list}
fastq-dump --outdir $GENESIS/ALIGN/454/FASTQ --defline-seq '@$sn[_$rn]/$ri' --gzip $GENESIS/ALIGN/454/SRA/${item}.sra
fastqc $GENESIS/ALIGN/454/FASTQ/${item}.fastq.gz --outdir $GENESIS/ALIGN/454/FASTQ/REPORTS --threads 20
trimmomatic SE -threads 20 -phred33 ${item}.fastq.gz ${item}_trimmed.fastq.gz HEADCROP:100 CROP 200 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 ILLUMINACLIP:$GENESIS/ALIGN/454/FASTQ/454_adapters.fa:2:30:10
fastqc $GENESIS/ALIGN/454/FASTQ/${item}_trimmed.fastq.gz --outdir $GENESIS/ALIGN/454/FASTQ/REPORTS2 --threads 20
done

gunzip *2.fastq.gz
find -type f -name '*2.fastq' -exec cat {} + > 454_all.fastq
gzip $GENESIS/ALIGN/454/FASTQ/454_all.fastq
fastqc $GENESIS/ALIGN/454/FASTQ/454_all.fastq.gz --outdir $GENESIS/ALIGN/454/FASTQ/REPORT_all --threads 40
Trinity --seqType fq --max_memory 192G --CPU 24 --single $GENESIS/ALIGN/454/FASTQ/454_all.fastq.gz --output $GENESIS/ALIGN/454/TRINITY_all/

## METASSEMBLY
Trinity --seqType fa --max_memory 192G --CPU 24 --single $GENESIS/ALIGN/454/TRINITY_all/Trinity.fasta,$GENESIS/ALIGN/MDV1/TRINITY/Trinity.fasta,$GENESIS/ALIGN/MDV23/TRINITY/Trinity.fasta,$GENESIS/ALIGN/VAD2/TRINITY/Trinity.fasta --output $GENESIS/ALIGN/TRINITY_FINAL_TRANSCRIPTOME/

# REDUCE BY USING PERL SCRIPT FROM TRINITY
perl /home/azken/Documentos/Miniconda3/envs/trinity/opt/trinity-2.8.4/util/misc/get_longest_isoform_seq_per_trinity_gene.pl Trinity2.fasta > ulmi_unigenes.fasta











