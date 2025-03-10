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

## CONTROL SAMPLES
/GENESIS/QUALITY/CTR$ fastqc /RAWDATA/CTR/*.fastq.gz --outdir REPORTS --threads 64
/GENESIS/QUALITY/CTR/REPORTS$ multiqc .

list='Gen-002 Gen-009 Gen-010 Gen-011 Gen-013 Gen-015 Gen-017 Gen-018 Gen-019 Gen-020 Gen-023 Gen-025 Gen-030 Gen-033 Gen-039 Gen-040 Gen-041 Gen-042 Gen-043 Gen-045 Gen-046 Gen-047 Gen-050 Gen-051 Gen-052 Gen-054 Gen-056 Gen-057 Gen-062 Gen-064 Gen-068 Gen-069 Gen-070 Gen-074 Gen-076 Gen-077 Gen-078 Gen-079 Gen-082 Gen-083 Gen-084 Gen-087 Gen-088 Gen-089 Gen-092 Gen-093 Gen-096 Gen-097 Gen-099 Gen-100 Gen-102 Gen-105 Gen-106 Gen-108 Gen-109 Gen-110 Gen-111 Gen-112 Gen-115 Gen-121 Gen-129 Gen-130 Gen-131 Gen-132 Gen-133 Gen-138 Gen-143 Gen-144 Gen-146 Gen-147 Gen-149 Gen-151 Gen-153 Gen-156 Gen-158 Gen-161 Gen-164 Gen-166 Gen-167 Gen-168 Gen-169 Gen-171 Gen-172 Gen-173 Gen-174 Gen-176 Gen-177 Gen-181 Gen-183 Gen-186 Gen-187 Gen-188 Gen-189 Gen-190 Gen-192 Gen-195'
for sample in ${list}
do
java -jar ~/Descargas/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 /RAWDATA/CTR/${sample}_1.fastq.gz /RAWDATA/CTR/${sample}_2.fastq.gz /TRIMMEDATA/CTR/${sample}_1_p.fastq.gz /TRIMMEDATA/CTR/${sample}_1_u.fastq.gz /TRIMMEDATA/CTR/${sample}_2_p.fastq.gz /TRIMMEDATA/CTR/${sample}_2_u.fastq.gz ILLUMINACLIP:TruSeq-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

/GENESIS/QUALITY/TRIMMED/CTR$ fastqc /TRIMMEDDATA/CTR/*.fastq.gz --outdir REPORTS --threads 64
/GENESIS/QUALITY/TRIMMED/CTR/REPORTS$ multiqc .

## INOCULATED SAMPLES
/GENESIS/QUALITY/EXP/$ fastqc /RAWDATA/EXP/*.fastq.gz --outdir REPORTS --threads 64
/GENESIS/QUALITY/EXP//REPORTS$ multiqc .

list='Gen-001 Gen-003 Gen-004 Gen-005 Gen-006 Gen-007 Gen-008 Gen-012 Gen-014 Gen-016 Gen-021 Gen-022 Gen-024 Gen-026 Gen-027 Gen-028 Gen-029 Gen-031 Gen-032 Gen-034 Gen-035 Gen-036 Gen-037 Gen-038 Gen-044 Gen-048 Gen-049 Gen-055 Gen-058 Gen-059 Gen-060 Gen-061 Gen-063 Gen-065 Gen-066 Gen-067 Gen-071 Gen-073 Gen-075 Gen-080 Gen-081 Gen-085 Gen-086 Gen-090 Gen-091 Gen-094 Gen-095 Gen-098 Gen-101 Gen-103 Gen-104 Gen-107 Gen-113 Gen-114 Gen-116 Gen-117 Gen-118 Gen-119 Gen-120 Gen-122 Gen-123 Gen-124 Gen-125 Gen-126 Gen-127 Gen-128 Gen-134 Gen-135 Gen-136 Gen-137 Gen-139 Gen-140 Gen-141 Gen-142 Gen-145 Gen-148 Gen-150 Gen-152 Gen-154 Gen-155 Gen-157 Gen-159 Gen-160 Gen-163 Gen-165 Gen-170 Gen-175 Gen-178 Gen-179 Gen-180 Gen-182 Gen-184 Gen-185 Gen-191 Gen-193 Gen-194'
for sample in ${list}
do
java -jar ~/Descargas/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 /RAWDATA/EXP/${sample}_1.fastq.gz /RAWDATA/EXP/${sample}_2.fastq.gz /TRIMMEDATA/EXP/${sample}_1_p.fastq.gz /TRIMMEDATA/EXP/${sample}_1_u.fastq.gz /TRIMMEDATA/EXP/${sample}_2_p.fastq.gz /TRIMMEDATA/EXP/${sample}_2_u.fastq.gz ILLUMINACLIP:TruSeq-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

/GENESIS/QUALITY/TRIMMED/EXP$ fastqc /TRIMMEDDATA/EXP/*.fastq.gz --outdir REPORTS --threads 64
/GENESIS/QUALITY/TRIMMED/EXP/REPORTS$ multiqc .

