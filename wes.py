#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Launcher script on a cluster HPC. 
from easyRun import *


with open('exome.cfg') as f:
    cfg = ConfigParser.ConfigParser()
    cfg.readfp(f)
QUEUE = cfg.get('sys', 'queue')
BINDIR = cfg.get('sys', 'bin')
REF = cfg.get('genome', 'ref')
BED = cfg.get('genome', 'bed')


# sys func
def createDir(path, mode='755', sudo=False):
    if not os.path.isdir(path):
        if sudo:
            os.system('sudo mkdir -m %s %s' % (mode, path))
        else:
            os.system('mkdir -m %s %s' % (mode, path))

def qsub(workDir, sample, cmd, afterID, queue=QUEUE, tag=None):  # tag: the same invoker invoke >1 qsub, etc. for different lanes, tag = LB 
    invoker = traceback.extract_stack()[-2][2]
    if tag != None:
        step = invoker + '_' + tag
    else:
        step = invoker
    shFile = '{workDir}/script/{sample}/{script}.sh'.format(workDir=workDir, sample=sample, script=step)
    step = sample + '\t' + step
    sh = script(cmd, afterID)
    logFile = '{workDir}/script/log'.format(workDir=workDir)
    return sh.qsub(shFile, queue, ppn=cfg.getint('ppn', invoker), logFile=logFile, logStep=step)

# gatk4.0 pip
def parseSampleLane(lane):
    tag = os.path.basename(lane)   # NP18E078_S5_L001
    m = re.match(r'(.+)_(S\d+_L\d+)', tag)
    ID, LB = m.group(1), m.group(2)
    return tag, ID, LB

def fqPreProc(workDir, sample, lane, afterID=None):  # lane=/lustre/rawData/analysis/wgs/180126_A00236_0031_AH732YDMXX/fq/NP18E078/NP18E078_S5_L001
    #support gz & gtz
    LB = parseSampleLane(lane)[2]
    cmds = []
    fq1, fq2 = sorted(glob.glob('%s_R*.fastq.gz*' % lane))
    if fq1.split('.')[-1] == 'gtz':
        cmds.append('gtz -c -d %s | gzip -c > %s.R1.fastq.gz' % (fq1, lane))
        cmds.append('rm %s' % fq1)
    else:
        cmds.append('mv %s %s.R1.fastq.gz' % (fq1, lane))
    if fq2.split('.')[-1] == 'gtz':
        cmds.append('gtz -c -d %s | gzip -c > %s.R2.fastq.gz' % (fq2, lane))
        cmds.append('rm %s' % fq2)
    else:
        cmds.append('mv %s %s.R2.fastq.gz' % (fq2, lane))
    cmds.append('md5sum {lane}.R1.fastq.gz > {lane}.R1.fastq.gz.md5'.format(lane=lane))
    cmds.append('md5sum {lane}.R2.fastq.gz > {lane}.R2.fastq.gz.md5'.format(lane=lane))
    cmds.append('fastqc {lane}.R1.fastq.gz'.format(lane=lane))
    cmds.append('fastqc {lane}.R2.fastq.gz'.format(lane=lane))
    cmd = '\n'.join(cmds)
    return qsub(workDir, sample, cmd, afterID, tag=LB)

def fqTrim(workDir, sample, lane, afterID=None):
    LB = parseSampleLane(lane)[2]
    cmds = []
    cmds.append('time java -jar /lustre/public/0.soft/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 20 -phred33 {lane}.R1.fastq.gz {lane}.R2.fastq.gz {lane}.R1.clean.fastq.gz {lane}.R1.unpair.fastq.gz {lane}.R2.clean.fastq.gz {lane}.R2.unpair.fastq.gz ILLUMINACLIP:/lustre/public/0.soft/Trimmomatic-0.36/adapters/TruSeq3-PE-3.fa:2:30:10:8:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:10 MINLEN:35'.format(lane=lane))
    cmds.append('rm {lane}.*unpair.fastq.gz'.format(lane=lane))
    cmd = '\n'.join(cmds)
    return qsub(workDir, sample, cmd, afterID, tag=LB)

def fqStat(workDir, sample, lane, afterID=None):
    LB = parseSampleLane(lane)[2]
    cmd = 'python {BINDIR}/fqStat.py {lane}'.format(BINDIR=BINDIR, lane=lane)
    return qsub(workDir, sample, cmd, afterID, tag=LB)

def bwa(workDir, sample, lane, outDir, afterID=None, ref=REF):
    tag, ID, LB = parseSampleLane(lane)
    cmd = 'time bwa mem -M -t 10 -R "@RG\\\\tID:{ID}\\\\tSM:{ID}\\\\tLB:{LB}\\\\tPL:ILLUMINA" {ref} {lane}.R1.clean.fastq.gz {lane}.R2.clean.fastq.gz | samtools view -S -b -o {outDir}/{tag}.raw.bam -'.format(ID=ID, LB=LB, lane=lane, outDir=outDir, tag=tag, ref=ref)
    return qsub(workDir, sample, cmd, afterID, tag=LB)

def sortBam(workDir, sample, lane, outDir, afterID=None, ref=REF):
    tag, _, LB = parseSampleLane(lane)
    cmds = []
    cmds.append('time gatk --java-options "-Xmx20g" SortSam --INPUT={outDir}/{tag}.raw.bam --OUTPUT={outDir}/{tag}.sort.bam --SORT_ORDER=coordinate --VALIDATION_STRINGENCY=SILENT'.format(outDir=outDir, tag=tag))
    cmds.append('time gatk --java-options "-Xmx20g" ReorderSam --INPUT={outDir}/{tag}.sort.bam --OUTPUT={outDir}/{tag}.reorder.bam --REFERENCE={ref} --VALIDATION_STRINGENCY=SILENT'.format(outDir=outDir, tag=tag, ref=ref))
    cmds.append('mv {outDir}/{tag}.reorder.bam {outDir}/{tag}.sort.bam'.format(outDir=outDir, tag=tag))
    cmd = '\n'.join(cmds)
    return qsub(workDir, sample, cmd, afterID, tag=LB)

def doLane(workDir, sample, lane, afterID=None):  # lane=/lustre/rawData/analysis/wgs/180126_A00236_0031_AH732YDMXX/fq/NP18E078/NP18E078_S5_L001
    outDir = '%s/alignment/%s' % (workDir, sample)
    job = fqPreProc(workDir, sample, lane, afterID)
    job = fqTrim(workDir, sample, lane, afterID=job)
    tmpjob = fqStat(workDir, sample, lane, afterID=job)
    job = bwa(workDir, sample, lane, outDir, afterID=job)
    job = sortBam(workDir, sample, lane, outDir, afterID=job)
    return [job, tmpjob]

def mergeLaneBams(workDir, sample, lanes=None, afterID=None):
    cmds = []
    if lanes:
        INPUT = ' '.join(map(lambda x: '--INPUT=%s/alignment/%s/%s.sort.bam' % (workDir, sample, x), lanes))
    else:
        INPUT = ' '.join(map(lambda x: '--INPUT=%s' % x, glob.glob('{workDir}/alignment/{sample}/*.sort.bam'.format(workDir=workDir, sample=sample))))
        #INPUT = glob.glob('{workDir}/alignment/{sample}/*.sort.bam'.format(workDir=workDir, sample=sample))
        #INPUT = map(lambda x: 'INPUT=%s' % x, INPUT)
        #INPUT = ' '.join(INPUT)
    OUTPUT = '--OUTPUT={workDir}/alignment/{sample}/{sample}.merge.bam'.format(workDir=workDir, sample=sample)
    cmds.append('time gatk --java-options "-Xmx30g" MergeSamFiles {INPUT} {OUTPUT}'.format(INPUT=INPUT, OUTPUT=OUTPUT))
    cmds.append('rm {workDir}/alignment/{sample}/*.raw.bam*'.format(workDir=workDir, sample=sample))
    cmds.append('rm {workDir}/alignment/{sample}/*.sort.bam*'.format(workDir=workDir, sample=sample))
    cmd = '\n'.join(cmds)
    return qsub(workDir, sample, cmd, afterID)

def dupmark(workDir, sample, afterID=None):
    cmds = []
    INPUT = '--INPUT={workDir}/alignment/{sample}/{sample}.merge.bam'.format(workDir=workDir, sample=sample)
    OUTPUT = '--OUTPUT={workDir}/alignment/{sample}/{sample}.dupmark.bam'.format(workDir=workDir, sample=sample)
    METRICS = '--METRICS_FILE={workDir}/alignment/{sample}/{sample}.dup.metrics'.format(workDir=workDir, sample=sample)
    cmds.append('time gatk --java-options "-Xmx20g" MarkDuplicates --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 {INPUT} {OUTPUT} {METRICS}'.format(INPUT=INPUT, OUTPUT=OUTPUT, METRICS=METRICS))
    cmds.append('rm {workDir}/alignment/{sample}/*.merge.bam*'.format(workDir=workDir, sample=sample))
    cmd = '\n'.join(cmds)
    return qsub(workDir, sample, cmd, afterID)

def bqsr(workDir, sample, afterID=None, ref=REF):
    cmds = []
    INPUT = '--input={workDir}/alignment/{sample}/{sample}.dupmark.bam'.format(workDir=workDir, sample=sample)
    OUTPUT = '--output={workDir}/alignment/{sample}/{sample}.recal_data.grp'.format(workDir=workDir, sample=sample)
    KNOWN = '--known-sites /lustre/public/1.DB/hg19/GATK_resource_bundle/dbsnp_138.hg19.vcf --known-sites /lustre/public/1.DB/hg19/GATK_resource_bundle/1000G_phase1.indels.hg19.sites.vcf --known-sites /lustre/public/1.DB/hg19/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'
    cmds.append('time gatk --java-options "-Xmx20g" BaseRecalibrator -R {ref} {INPUT} {KNOWN} {OUTPUT}'.format(ref=ref, KNOWN=KNOWN, INPUT=INPUT, OUTPUT=OUTPUT))
    BQSR = '--bqsr-recal-file={workDir}/alignment/{sample}/{sample}.recal_data.grp'.format(workDir=workDir, sample=sample)
    OUTPUT = '--output={workDir}/alignment/{sample}/{sample}.final.bam'.format(workDir=workDir, sample=sample)
    cmds.append('time gatk --java-options "-Xmx20g" ApplyBQSR -R {ref} {INPUT} {BQSR} {OUTPUT}'.format(ref=ref, BQSR=BQSR, INPUT=INPUT, OUTPUT=OUTPUT))
    cmds.append('rm {workDir}/alignment/{sample}/*.dupmark.bam*'.format(workDir=workDir, sample=sample))
    cmd = '\n'.join(cmds)
    return qsub(workDir, sample, cmd, afterID)

def haplotype(workDir, sample, afterID=None, ref=REF, bed=BED):
    INPUT = '--input={workDir}/alignment/{sample}/{sample}.final.bam'.format(workDir=workDir, sample=sample)
    OUTPUT = '--output={workDir}/variation/{sample}/{sample}.g.vcf.gz'.format(workDir=workDir, sample=sample)
    cmd = 'time gatk --java-options "-Xmx20g" HaplotypeCaller -R {ref} --emit-ref-confidence GVCF -L {bed} {INPUT} {OUTPUT}'.format(ref=ref, bed=bed, INPUT=INPUT, OUTPUT=OUTPUT)
    return qsub(workDir, sample, cmd, afterID)

def combineGVCFs(workDir, samples, afterID=None, ref=REF, L=50):
    createDir('%s/script/cohort' % workDir)
    createDir('%s/variation/cohort' % workDir)
    VAR = ' '.join(map(lambda x: '--variant=%s/variation/%s/%s.g.vcf.gz' % (workDir, x, x), samples))
    if len(samples) < L:
        N = L - len(samples)
        VAR = ' '.join([VAR, add_gvcf(N)])
    OUTPUT = '--output={workDir}/variation/cohort/cohort.g.vcf.gz'.format(workDir=workDir)
    cmd = 'time gatk --java-options "-Xmx100g" CombineGVCFs -R {ref} {VAR} {OUTPUT}'.format(ref=ref, VAR=VAR, OUTPUT=OUTPUT)
    return qsub(workDir, 'cohort', cmd, afterID)

def genotype(workDir, sample, afterID=None, ref=REF):
    VAR = '--variant={workDir}/variation/{sample}/{sample}.g.vcf.gz'.format(workDir=workDir, sample=sample)
    OUTPUT = '--output={workDir}/variation/{sample}/{sample}.raw.vcf'.format(workDir=workDir, sample=sample)
    DBSNP = '--dbsnp=/lustre/public/1.DB/hg19/GATK_resource_bundle/dbsnp_138.hg19.vcf'
    cmd = 'time gatk --java-options "-Xmx20g" GenotypeGVCFs -R {ref} {VAR} {DBSNP} {OUTPUT}'.format(ref=ref, VAR=VAR, DBSNP=DBSNP, OUTPUT=OUTPUT)
    return qsub(workDir, sample, cmd, afterID)

def vqsr_snp(workDir, sample, afterID=None, ref=REF, InbreedingCoeff=True):
    cmds = []
    cmds.append('time gatk SelectVariants -R {ref} -V {workDir}/variation/{sample}/{sample}.raw.vcf --select-type-to-include SNP -O {workDir}/variation/{sample}/{sample}.snp.vcf'.format(ref=ref, workDir=workDir, sample=sample))
    VAR = '--variant={workDir}/variation/{sample}/{sample}.snp.vcf'.format(workDir=workDir, sample=sample)
    RES = '--resource hapmap,known=false,training=true,truth=true,prior=15.0:/lustre/public/1.DB/hg19/GATK_resource_bundle/hapmap_3.3.hg19.sites.vcf --resource omni,known=false,training=true,truth=false,prior=12.0:/lustre/public/1.DB/hg19/GATK_resource_bundle/1000G_omni2.5.hg19.sites.vcf --resource 1000G,known=false,training=true,truth=false,prior=10.0:/lustre/public/1.DB/hg19/GATK_resource_bundle/1000G_phase1.snps.high_confidence.hg19.sites.vcf --resource dbsnp,known=true,training=false,truth=false,prior=2.0:/lustre/public/1.DB/hg19/GATK_resource_bundle/dbsnp_138.hg19.vcf'
    if InbreedingCoeff:
        ANNO = ' '.join(map(lambda x: '-an %s' % x, ['QD', 'MQ', 'MQRankSum', 'ReadPosRankSum', 'FS', 'SOR', 'InbreedingCoeff']))  # InbreedingCoeff needs at least 10 samples
    else:
        ANNO = ' '.join(map(lambda x: '-an %s' % x, ['QD', 'MQ', 'MQRankSum', 'ReadPosRankSum', 'FS', 'SOR']))
    TRANCH = '--tranches-file={workDir}/variation/{sample}/{sample}.snp.tranches'.format(workDir=workDir, sample=sample)
    OUTPUT = '--output={workDir}/variation/{sample}/{sample}.snp.recal'.format(workDir=workDir, sample=sample)
    cmds.append('time gatk --java-options "-Xmx20g" VariantRecalibrator -R {ref} -mode SNP {VAR} {RES} {ANNO} --max-gaussians 4 {TRANCH} {OUTPUT}'.format(ref=ref, VAR=VAR, RES=RES, ANNO=ANNO, TRANCH=TRANCH, OUTPUT=OUTPUT))
    RECAL = '--recal-file={workDir}/variation/{sample}/{sample}.snp.recal'.format(workDir=workDir, sample=sample)
    OUTPUT = '--output={workDir}/variation/{sample}/{sample}.snp_recal.vcf'.format(workDir=workDir, sample=sample)
    cmds.append('time gatk --java-options "-Xmx20g" ApplyVQSR -R {ref} -mode SNP {VAR} {RECAL} -ts-filter-level 99.0 {TRANCH} {OUTPUT}'.format(ref=ref, VAR=VAR, RECAL=RECAL, TRANCH=TRANCH, OUTPUT=OUTPUT))
    cmd = '\n'.join(cmds)
    return qsub(workDir, sample, cmd, afterID)

def vqsr_indel(workDir, sample, afterID=None, ref=REF, InbreedingCoeff=True):
    cmds = []
    cmds.append('time gatk SelectVariants -R {ref} -V {workDir}/variation/{sample}/{sample}.raw.vcf --select-type-to-include INDEL -O {workDir}/variation/{sample}/{sample}.indel.vcf'.format(ref=ref, workDir=workDir, sample=sample))
    VAR = '--variant={workDir}/variation/{sample}/{sample}.indel.vcf'.format(workDir=workDir, sample=sample)
    RES = '--resource mills,known=false,training=true,truth=true,prior=12.0:/lustre/public/1.DB/hg19/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf --resource dbsnp,known=true,training=false,truth=false,prior=2.0:/lustre/public/1.DB/hg19/GATK_resource_bundle/dbsnp_138.hg19.vcf'
    if InbreedingCoeff:
        ANNO = ' '.join(map(lambda x: '-an %s' % x, ['QD', 'DP', 'FS', 'SOR', 'ReadPosRankSum', 'MQRankSum', 'InbreedingCoeff']))  #note that 'DP' replace 'MQ' rather than SNP
    else:
        ANNO = ' '.join(map(lambda x: '-an %s' % x, ['QD', 'DP', 'FS', 'SOR', 'ReadPosRankSum', 'MQRankSum']))
    TRANCH = '--tranches-file={workDir}/variation/{sample}/{sample}.indel.tranches'.format(workDir=workDir, sample=sample)
    OUTPUT = '--output={workDir}/variation/{sample}/{sample}.indel.recal'.format(workDir=workDir, sample=sample)
    cmds.append('time gatk --java-options "-Xmx20g" VariantRecalibrator -R {ref} -mode INDEL {VAR} {RES} {ANNO} --max-gaussians 4 {TRANCH} {OUTPUT}'.format(ref=ref, VAR=VAR, RES=RES, ANNO=ANNO, TRANCH=TRANCH, OUTPUT=OUTPUT))
    RECAL = '--recal-file={workDir}/variation/{sample}/{sample}.indel.recal'.format(workDir=workDir, sample=sample)
    OUTPUT = '--output={workDir}/variation/{sample}/{sample}.indel_recal.vcf'.format(workDir=workDir, sample=sample)
    cmds.append('time gatk --java-options "-Xmx20g" ApplyVQSR -R {ref} -mode INDEL {VAR} {RECAL} -ts-filter-level 99.0 {TRANCH} {OUTPUT}'.format(ref=ref, VAR=VAR, RECAL=RECAL, TRANCH=TRANCH, OUTPUT=OUTPUT))
    cmd = '\n'.join(cmds)
    return qsub(workDir, sample, cmd, afterID)

def vqsr(workDir, sample, afterID=None, ref=REF):
    job = []
    job.append(vqsr_snp(workDir, sample, afterID))
    job.append(vqsr_indel(workDir, sample, afterID))
    cmds = []
    cmds.append('time gatk SelectVariants -R {ref} -V {workDir}/variation/{sample}/{sample}.raw.vcf --select-type-to-include MIXED -O {workDir}/variation/{sample}/{sample}.mix.vcf'.format(ref=ref, workDir=workDir, sample=sample))
    cmds.append('time gatk MergeVcfs -I {workDir}/variation/{sample}/{sample}.snp_recal.vcf -I {workDir}/variation/{sample}/{sample}.indel_recal.vcf -I {workDir}/variation/{sample}/{sample}.mix.vcf -O {workDir}/variation/{sample}/{sample}.final.vcf'.format(workDir=workDir, sample=sample))
    cmds.append('rm {workDir}/variation/{sample}/*snp*'.format(workDir=workDir, sample=sample))
    cmds.append('rm {workDir}/variation/{sample}/*indel*'.format(workDir=workDir, sample=sample))
    cmds.append('rm {workDir}/variation/{sample}/*mix*'.format(workDir=workDir, sample=sample))
    cmds.append('rm {workDir}/variation/{sample}/*raw*'.format(workDir=workDir, sample=sample))
    cmd = '\n'.join(cmds)
    return qsub(workDir, sample, cmd, job)

def hardFilter(workDir, sample, afterID=None, ref=REF):
    cmds = []
    cmds.append('time gatk SelectVariants -R {ref} -V {workDir}/variation/{sample}/{sample}.raw.vcf --select-type-to-include SNP -O {workDir}/variation/{sample}/{sample}.snp.vcf'.format(ref=ref, workDir=workDir, sample=sample))
    VAR = '--variant={workDir}/variation/{sample}/{sample}.snp.vcf'.format(workDir=workDir, sample=sample)
    OUTPUT = '--output={workDir}/variation/{sample}/{sample}.snp.filter.vcf'.format(workDir=workDir, sample=sample)
    cmds.append('time gatk --java-options "-Xmx20g" VariantFiltration -R {ref} {VAR} {OUTPUT} --filter-expression \"DP<=0 || QD<2.0 || MQ<40.0 || FS>60.0 || HaplotypeScore>60.0 ||  MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --filter-name \"StandardFilterSNV\"'.format(ref=ref, VAR=VAR, OUTPUT=OUTPUT))
    
    cmds.append('time gatk SelectVariants -R {ref} -V {workDir}/variation/{sample}/{sample}.raw.vcf --select-type-to-include INDEL -O {workDir}/variation/{sample}/{sample}.indel.vcf'.format(ref=ref, workDir=workDir, sample=sample))
    VAR = '--variant={workDir}/variation/{sample}/{sample}.indel.vcf'.format(workDir=workDir, sample=sample)
    OUTPUT = '--output={workDir}/variation/{sample}/{sample}.indel.filter.vcf'.format(workDir=workDir, sample=sample)
    cmds.append('time gatk --java-options "-Xmx20g" VariantFiltration -R {ref} {VAR} {OUTPUT} --filter-expression \"DP<=0 || FS>200.0 || ReadPosRankSum < -20.0\" --filter-name \"StandardFilterINDEL\"'.format(ref=ref, VAR=VAR, OUTPUT=OUTPUT))

    cmds.append('time gatk SelectVariants -R {ref} -V {workDir}/variation/{sample}/{sample}.raw.vcf --select-type-to-include MIXED -O {workDir}/variation/{sample}/{sample}.mix.vcf'.format(ref=ref, workDir=workDir, sample=sample))
    cmds.append('time gatk MergeVcfs -I {workDir}/variation/{sample}/{sample}.snp.filter.vcf -I {workDir}/variation/{sample}/{sample}.indel.filter.vcf -I {workDir}/variation/{sample}/{sample}.mix.vcf -O {workDir}/variation/{sample}/{sample}.final.vcf'.format(workDir=workDir, sample=sample))

    cmds.append('rm {workDir}/variation/{sample}/*snp*'.format(workDir=workDir, sample=sample))
    cmds.append('rm {workDir}/variation/{sample}/*indel*'.format(workDir=workDir, sample=sample))
    cmds.append('rm {workDir}/variation/{sample}/*mix*'.format(workDir=workDir, sample=sample))
    cmds.append('rm {workDir}/variation/{sample}/*raw*'.format(workDir=workDir, sample=sample))
    cmd = '\n'.join(cmds)
    return qsub(workDir, sample, cmd, afterID)

################ add-on #####################
def std_id(sample):
    m = re.search(r'(t\w+)-(\w+)', sample)
    if m:
        return m.group(1) + '.' + m.group(2)
    else:
        return 'tEXOME2.' + sample

def std_ids(samples):
    IDs = []
    for sample in samples:
        IDs.append(std_id(sample))
    return IDs

def depthQC(workDir, sample, afterID=None, bed='/lustre/public/1.DB/panel/tEXOME2.bed_extend_10bp'):
    with open('%s/script/%s/depthQC.sh' % (workDir, sample), 'w') as out:
        out.write('perl /lustre/usr/bin/depthQC/depth.pl -m {workDir}/alignment/{sample}/{sample}.final.bam -d /lustre/public/1.DB/hg19/hg19.dict -o {workDir}/depth_QC/{sample} -b {bed}\n'.format(workDir=workDir, sample=sample, bed=bed))
        out.write(log(workDir, sample, 'depthQC'))
        out.write('Rscript /lustre/usr/bin/depthQC/dis_low.R {workDir}/depth_QC/{sample}.perBase.depth {sample} {workDir}/depth_QC/{sample}.perBase.png\n'.format(workDir=workDir, sample=sample))
        out.write('Rscript /lustre/usr/bin/depthQC/cumn.R {workDir}/depth_QC/{sample}.Cumu.depth {sample} {workDir}/depth_QC/{sample}.Cumu.png\n'.format(workDir=workDir, sample=sample))
        out.write('Rscript /lustre/usr/bin/depthQC/insertsize.R {workDir}/depth_QC/{sample}.insertsize {sample} {workDir}/depth_QC/{sample}.insertsize.png\n'.format(workDir=workDir, sample=sample))
    return qsub('{workDir}/script/{sample}/depthQC.sh'.format(workDir=workDir, sample=sample), afterID)

def coverCnt(workDir, sample, afterID=None, bed=BED):
    ID = std_id(sample)
    with open('%s/script/%s/coverCnt.sh' % (workDir, sample), 'w') as out:
        out.write('/lustre/public/0.soft/bedtools-2.17.0/bin/coverageBed -abam {workDir}/alignment/{sample}/{sample}.final.bam -hist -d -b {bed} > {workDir}/annotation/coverage/{ID}.coverage\n'.format(workDir=workDir, sample=sample, bed=bed, ID=ID))
        out.write('python /lustre/usr/bin/coverCnt.py {workDir}/annotation/coverage/{ID}.coverage\n'.format(workDir=workDir, ID=ID))
        out.write(log(workDir, sample, 'coverCnt'))
        out.write('gzip {workDir}/annotation/coverage/{ID}.coverage\n'.format(workDir=workDir, ID=ID))
    return qsub('{workDir}/script/{sample}/coverCnt.sh'.format(workDir=workDir, sample=sample), afterID)
#############################################

def add_gvcf(k):
    gvcfs = []
    with open('/lustre/rawData/analysis/zzz/exome/add_gvcf.list') as f:
        for line in f.readlines():
            gvcfs.append(line.strip())
    return ' '.join(map(lambda x: '--variant=%s' % x, gvcfs[:k]))

def splitCohortVcf(workDir, sample, selectSamples, afterID=None):
    with open('%s/script/%s/splitCohortVcf.sh' % (workDir, sample), 'w') as out:
        if selectSamples == None:
            out.write('python /lustre/public/test/liujx/wes/splitCohortVcf.py -vcf {workDir}/variation/{sample}/{sample}.final.vcf\n'.format(workDir=workDir, sample=sample))
        else:
            out.write('python /lustre/public/test/liujx/wes/splitCohortVcf.py -vcf {workDir}/variation/{sample}/{sample}.final.vcf -select {selectSamples}\n'.format(workDir=workDir, sample=sample, selectSamples=','.join(selectSamples)))
        out.write(log(workDir, sample, 'splitCohortVcf'))
    return qsub('{workDir}/script/{sample}/splitCohortVcf.sh'.format(workDir=workDir, sample=sample), afterID)

def runSample(workDir, sample, zipType='gz'):
    createDir('%s/fq/%s' % (workDir, sample))
    createDir('%s/script/%s' % (workDir, sample))
    createDir('%s/alignment/%s' % (workDir, sample))
    createDir('%s/variation/%s' % (workDir, sample))
    jobid = []
    lanes = []
    if zipType == 'gtz':
        fqs = sorted(glob.glob('{workDir}/bcl2fastq/{sample}/*_R1_*.fastq.gz.gtz'.format(workDir=workDir, sample=sample)))
    elif zipType == 'gz':
        fqs = sorted(glob.glob('{workDir}/bcl2fastq/{sample}/*_R1_*.fastq.gz'.format(workDir=workDir, sample=sample)))
    for fq in fqs:
        m = re.match(r'(.+)_R1_.*.fastq.gz', fq)
        lane = m.group(1)
        jobid.append(doLane(workDir, sample, lane, zipType))  # lane=/lustre/rawData/analysis/zzz/exome/170331_NB501851_0045_AH7LYTBGX2/fq/KY18K1279/KY18K1279_S1_L001
        lane = re.sub(r'bcl2fastq', 'fq', lane)
        lanes.append(os.path.basename(lane))  # lane=/lustre/rawData/analysis/zzz/exome/170331_NB501851_0045_AH7LYTBGX2/bcl2fastq/KY18K1279/KY18K1279_S1_L001
    jobid = merge(workDir, sample, lanes, jobid)
    jobid = dupmark(workDir, sample, jobid)
    jobid = bqsr(workDir, sample, jobid)
    # final bam, get!
    if sample[:3] != 'NTC':
        qc1 = depthQC(workDir, sample, jobid)
        qc2 = coverCnt(workDir, sample, jobid)
    # run single sample
    jobid = haplotype(workDir, sample, jobid)
    jobid = genotype(workDir, sample, jobid)
    jobid = hardFilter(workDir, sample, jobid)  # jobid = vqsr(workDir, sample, jobid, InbreedingCoeff=False)
    jobid = splitCohortVcf(workDir, sample, None, jobid)
    return [jobid, qc1, qc2] if sample[:3] != 'NTC' else jobid

def select_cat_snp_indel(workDir, samples, afterID=None):
    with open('%s/script/select_cat_snp_indel.sh' % workDir, 'w') as out:
        out.write('python /lustre/public/test/liujx/wes/select_cat_snp_indel.py {workDir} {samples}\n'.format(workDir=workDir, samples=','.join(samples)))
        out.write(log(workDir, 'cohort', 'select_cat_snp_indel'))
    return qsub('{workDir}/script/select_cat_snp_indel.sh'.format(workDir=workDir), afterID)

def runCohort(workDir, samples, afterID=None):
    afterID = combineGVCFs(workDir, samples, afterID)
    afterID = genotype(workDir, 'cohort', afterID)
    afterID = vqsr(workDir, 'cohort', afterID)
    afterID = splitCohortVcf(workDir, 'cohort', samples, afterID)
    afterID = select_cat_snp_indel(workDir, samples, afterID)
    return afterID


# cnv pip, including chr ploidy, etc 21T,XXY
def exonCNV(workDir, afterID=None, config=config['exonCNV'], bed='/lustre/public/1.DB/panel/tEXOME2.bed'):
    with open('%s/script/exonCNV.sh' % workDir, 'w') as out:
        out.write('perl /lustre/home/sbsuser/script/exonCNV/run_exonCNV.pl {workDir} {workDir}/exonCNV/ {bed}\n'.format(workDir=workDir, bed=bed))
        out.write('perl /lustre/usr/bin/genderChr.pl {workDir}\n'.format(workDir=workDir))
        out.write('chmod 750 {workDir}/exonCNV/script\n'.format(workDir=workDir))
        out.write(log(workDir, 'cohort', 'exonCNV'))
    return qsub('{workDir}/script/exonCNV.sh'.format(workDir=workDir), afterID, config)


# annotation pip
def std_merge_vcf(workDir, annoDir, samples, IDs, afterID=None):
    with open('%s/script/std_merge_vcf.sh' % workDir, 'w') as out:
        out.write('python /lustre/public/test/liujx/wes/std_merge_vcf.py {workDir} {annoDir} {samples} {IDs}\n'.format(workDir=workDir, annoDir=annoDir, samples=','.join(samples), IDs=','.join(IDs)))
        out.write(log(workDir, 'cohort', 'std_merge_vcf'))
    return qsub('{workDir}/script/std_merge_vcf.sh'.format(workDir=workDir), afterID)

def vcf2csv(workDir, annoDir, sample, ID, afterID=None, platform='exome', config=config['vcf2csv']):
    vcf = '{annoDir}/input/{ID}.vcf.convert'.format(annoDir=annoDir, ID=ID)
    ID = ID.split('.')[-1]
    with open('%s/script/%s/vcf2csv.sh' % (workDir, sample), 'w') as out:
        out.write('time python {binDir}/vcf2annovar.pyc -infile {vcf} -outPrefix {annoDir}/rawAnnovar/{ID} &>> {annoDir}/logs/{ID}.log\n'.format(binDir=binDir, vcf=vcf, annoDir=annoDir, ID=ID))
        out.write(log(workDir, sample, 'vcf2annovar'))
        createDir('%s/%s' % (annoDir, ID))
        out.write('time python {binDir}/varClassify.pyc -infile {annoDir}/rawAnnovar/{ID}.genome_summary.csv -platform {platform} -batchfile {annoDir}/merge.tmp -outfile {annoDir}/{ID}/{ID}.csv &>> {annoDir}/logs/{ID}.log\n'.format(binDir=binDir, annoDir=annoDir, ID=ID, platform=platform))
        out.write(log(workDir, sample, 'varClassify'))
        out.write('time python {binDir}/varStat.py {annoDir}/{ID}/{ID}.csv &>> {annoDir}/logs/{ID}.log\n'.format(binDir=binDir, annoDir=annoDir, ID=ID))
        #out.write(log(workDir, sample, 'varStat'))
        out.write('chmod 750 {annoDir}/{ID}/{ID}.csv\n'.format(ID=ID, annoDir=annoDir))
        out.write('perl {binDir}/varAdjacent.pl ^_^ {ID} {annoDir}\n'.format(binDir=binDir, ID=ID, annoDir=annoDir))
    return qsub('{workDir}/script/{sample}/vcf2csv.sh'.format(workDir=workDir, sample=sample), afterID, config)

def annotation(workDir, samples, IDs, afterID=None):
    annoDir = '%s/annotation' % workDir
    createDir(annoDir)
    createDir('%s/input' % annoDir, mode='754')
    createDir('%s/rawAnnovar' % annoDir, mode='750')
    createDir('%s/logs' % annoDir, mode='750')
    afterID = std_merge_vcf(workDir, annoDir, samples, IDs, afterID)
    jobid = []
    for sample,ID in zip(samples, IDs):
        jobid.append(vcf2csv(workDir, annoDir, sample, ID, afterID))
    return jobid


# qual ctrl
def coverStat(workDir, afterID=None):  # {workDir}/annotation/coverage中glob coverCnt产生的 stat文件
    with open('%s/script/coverStat.sh' % workDir, 'w') as out:
        out.write('python /lustre/usr/bin/coverStat.py {workDir}/annotation/coverage {batch}\n'.format(workDir=workDir, batch=os.path.basename(workDir)))
        out.write(log(workDir, 'cohort', 'coverStat'))
    return qsub('{workDir}/script/coverStat.sh'.format(workDir=workDir), afterID)

def qcReport(workDir, afterID=None):
    with open('%s/script/qcReport.sh' % workDir, 'w') as out:
        out.write('python /lustre/usr/bin/qcReport.py %s\n' % workDir)
        out.write(log(workDir, 'cohort', 'qcReport'))
    return qsub('{workDir}/script/qcReport.sh'.format(workDir=workDir), afterID)

def mergeXls(workDir, sampleCnt, afterID=None):
    ppn = 8 if sampleCnt > 100 else 4
    with open('%s/script/mergeXls.sh' % workDir, 'w') as out:
        out.write('python /lustre/usr/bin/mergeXls.py %s\n' % workDir)
        out.write(log(workDir, 'cohort', 'mergeXls'))
    return qsub('{workDir}/script/mergeXls.sh'.format(workDir=workDir), afterID, config='-l nodes=1:ppn=%d' % ppn)



if __name__ == "__main__":
    main(sys.argv[1:])