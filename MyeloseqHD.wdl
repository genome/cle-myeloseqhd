version 1.0

import "MyeloseqHDAnalysis.wdl" as subWF

workflow MyeloseqHD {

    input {

        File SampleSheet
        # sample sheet has this structure:
        # index  name  RG_ID  RG_FLOWCELL  RG_LANE  RG_LIB  RG_SAMPLE MRN ALL_MRNs ACCESSION DOB SEX EXCEPTION [R1] [R2]

        File? DemuxSampleSheet
        String? IlluminaDir
        String? DragenEnv

        String JobGroup
        String OutputDir
        String RunInfoString

        Boolean DataTransfer
        String? CasesExcluded

        String Queue
        String DragenQueue
        String DragenDockerImage
        String VariantDB

        Int MinReads
        Float MinVaf
        Int readfamilysize = 3

        Array[String] Adapters = ["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC","AGATCGGAAGAGCGTCGTGTAGGGAAA"]

        String DragenReference = "/storage1/fs1/gtac-mgi/Active/CLE/reference/dragen424_hg38"
        String Reference       = "/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/refdata/hg38/all_sequences.fa"
        String ReferenceDict   = "/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/refdata/hg38/all_sequences.dict"
        String VEP             = "/storage1/fs1/gtac-mgi/Active/CLE/reference/VEP_cache"

        String MyeloSeqHDRepo
    }

    String CoverageBed = MyeloSeqHDRepo + "/accessory_files/MyeloseqHD.CoverageQC.hg38.bed"
    String Hotspot     = MyeloSeqHDRepo + "/accessory_files/MyeloseqHD.hotspots.vcf.gz"

    String QC_pl   = MyeloSeqHDRepo + "/dockerfiles/docker-myeloseq/QC_metrics.pl"
    String xfer_pl = MyeloSeqHDRepo + "/scripts/data_transfer.pl"

    String DemuxFastqDir = "/storage1/fs1/gtac-mgi/Active/CLE/assay/myeloseqhd/demux_fastq"


    if (defined(DemuxSampleSheet)){
      call dragen_demux {
        input: Dir=IlluminaDir,
        OutputDir=OutputDir,
        SampleSheet=DemuxSampleSheet,
        DragenEnv=DragenEnv,
        DragenDockerImage=DragenDockerImage,
        queue=DragenQueue,
        jobGroup=JobGroup
      }

      call prepare_samples {
        input: SampleSheet=SampleSheet,
        Fastq1=dragen_demux.read1,
        Fastq2=dragen_demux.read2,
        queue=Queue,
        jobGroup=JobGroup
      }
    }

    Array[Array[String]] inputData = read_tsv(select_first([prepare_samples.sample_sheet,SampleSheet]))

    # the inputdata should be: index  name  RG_ID  RG_FLOWCELL  RG_LANE  RG_LIB  RG_SAMPLE MRN ALL_MRNs ACCESSION DOB SEX EXCEPTION read1path read2path
    scatter (samples in inputData){

        if(!defined(DemuxSampleSheet)){
          call trim_adapters {
              input: Read1=samples[13],
              Read2=samples[14],
              Adapters=Adapters,
              Name=samples[1],
              queue=Queue,
              jobGroup=JobGroup
          }
        }

        call dragen_align {
            input: DragenRef=DragenReference,
                   Hotspot=Hotspot,
                   fastq1=select_first([trim_adapters.read1,samples[13]]),
                   fastq2=select_first([trim_adapters.read2,samples[14]]),
                   Name=samples[1],
                   RG=samples[3] + '.' + samples[4] + '.' + samples[0],
                   SM=samples[6],
                   LB=samples[5] + '.' + samples[0],
                   readfamilysize=readfamilysize,
                   CoverageBed=CoverageBed,
                   OutputDir=OutputDir,
                   SubDir=samples[1] + '_' + samples[0],
                   DragenEnv=DragenEnv,
                   DragenDockerImage=DragenDockerImage,
                   queue=DragenQueue,
                   jobGroup=JobGroup
        }

        call subWF.MyeloseqHDAnalysis {
            input: Bam=dragen_align.bam,
                   BamIndex=dragen_align.bai,
                   DragenVcf=dragen_align.vcf,
                   DragenVcfIndex=dragen_align.index,
                   refFasta=Reference,
                   ReferenceDict=ReferenceDict,
                   Name=samples[1],
                   mrn=samples[7],
                   all_mrn=samples[8],
                   accession=samples[9],
                   DOB=samples[10],
                   sex=samples[11],
                   exception=samples[12],
                   RunInfoString=RunInfoString,
                   VariantDB=VariantDB,
                   Vepcache=VEP,
                   MinReads=MinReads,
                   MinVaf=MinVaf,
                   MyeloSeqHDRepo=MyeloSeqHDRepo,
                   OutputDir=OutputDir,
                   SubDir=samples[1] + '_' + samples[0],
                   Queue=Queue,
                   JobGroup=JobGroup
        }
    }

    if (defined(DemuxSampleSheet)){
        call move_demux_fastq {
            input: order_by=MyeloseqHDAnalysis.all_done,
            Batch=basename(OutputDir),
            DemuxFastqDir=DemuxFastqDir,
            queue=DragenQueue,
            jobGroup=JobGroup
        }
    }

    call batch_qc {
        input: order_by=MyeloseqHDAnalysis.all_done,
               BatchDir=OutputDir,
               QC_pl=QC_pl,
               queue=Queue,
               jobGroup=JobGroup
    }

    if (DataTransfer) {
        call data_transfer {
            input: order_by=batch_qc.done,
            BatchDir=OutputDir,
            xfer_pl=xfer_pl,
            excluded=CasesExcluded,
            queue=Queue,
            jobGroup=JobGroup
        }
    }
}


task dragen_demux {
     input {
         String? Dir
         File? SampleSheet
         String? DragenEnv
         String OutputDir
         String jobGroup
         String queue
         String DragenDockerImage
     }

     String batch = basename(OutputDir)
     String StagingDir = "/staging/runs/MyeloSeqHD/"
     String LocalFastqDir = StagingDir + "demux_fastq/" + batch
     String LocalReportDir = LocalFastqDir + "/Reports"
     String LocalSampleSheet = StagingDir + "sample_sheet/" + batch + '.csv'
     String log = StagingDir + "log/" + batch + "_demux.log"
     String DemuxReportDir = OutputDir + "/dragen_demux_reports"

     command {
         /bin/cp ${SampleSheet} ${LocalSampleSheet} && \
         /opt/edico/bin/dragen --bcl-conversion-only true --bcl-only-matched-reads true --strict-mode true --sample-sheet ${LocalSampleSheet} --bcl-input-directory ${Dir} --output-directory ${LocalFastqDir} &> ${log} && \
         /bin/ls ${LocalFastqDir}/*_R1_001.fastq.gz > Read1_list.txt && \
         /bin/ls ${LocalFastqDir}/*_R2_001.fastq.gz > Read2_list.txt && \
         /bin/mv ${log} ./ && \
         /bin/rm -f ${LocalSampleSheet} && \
         /bin/cp -r ${LocalReportDir} ${DemuxReportDir} && \

         if [ -n "$(/bin/find ${LocalFastqDir}/TW*.fastq.gz -type f -size -10M)" ]; then
             echo 'demux fastq size checking failed !' && \
             exit 1
         fi
     }

     runtime {
         docker_image: DragenDockerImage
         dragen_env: DragenEnv
         cpu: "20"
         memory: "200 G"
         queue: queue
         job_group: jobGroup 
     }
     output {
         File read1 = "Read1_list.txt"
         File read2 = "Read2_list.txt"
     }
}

task prepare_samples {
     input {
         File SampleSheet
         String Fastq1
         String Fastq2
         String jobGroup
         String queue
     }

     command <<<
             /bin/cp ~{Fastq1} 1.tmp.txt
             /bin/cp ~{Fastq2} 2.tmp.txt
             /usr/bin/perl -e 'open(R1,"1.tmp.txt"); @r1 = <R1>; \
                 chomp @r1; close R1;\
                 open(R2,"2.tmp.txt"); @r2 = <R2>; \
                 chomp @r2; close R2; \
                 open(SS,"~{SampleSheet}");
                 while(<SS>){
                     chomp;
                     my @l = split("\t",$_);
                     my $s = $l[1].'_';
                     my $r1 = (grep /$s/, @r1)[0];
                     my $r2 = (grep /$s/, @r2)[0];
                     print join("\t",@l,$r1,$r2),"\n";
                 }
                 close SS;' > sample_sheet.txt
     >>>
     runtime {
         docker_image: "docker1(ubuntu:xenial)"
         cpu: "1"
         memory: "4 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File sample_sheet = "sample_sheet.txt"
     }
}

task trim_adapters {
     input{
         String Read1
         String Read2
         Array[String] Adapters
         String Name
         String jobGroup
         String queue
     }
     command {
        if [[ "${Read1}" == *"_R1_001"* ]]; then
             /bin/cp ${Read1} ${Name}.1.fastq.gz && \
             /bin/cp ${Read2} ${Name}.2.fastq.gz
        else
             export PYTHONPATH=/opt/cutadapt/lib/python2.7/site-packages/ && \
             /opt/cutadapt/bin/cutadapt -a ${Adapters[0]} -A ${Adapters[1]} -o ${Name}.1.fastq.gz -p ${Name}.2.fastq.gz ${Read1} ${Read2}
        fi
     }

     runtime {
         docker_image: "docker1(dhspence/docker-cutadapt:latest)"
         cpu: "1"
         memory: "8 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File read1 = "${Name}.1.fastq.gz"
         File read2 = "${Name}.2.fastq.gz"
     }
}

task dragen_align {
     input{
         String Name
         String DragenRef
         String Hotspot
         String fastq1
         String fastq2
         String RG
         String SM
         String LB
         String CoverageBed
         String OutputDir
         String SubDir
         String jobGroup
         String queue
         String DragenDockerImage
         String? DragenEnv

         Int? TrimLen
         Int readfamilysize
     }

     String batch = basename(OutputDir)
     String StagingDir = "/staging/runs/MyeloSeqHD/"
     String LocalAlignDir = StagingDir + "align/" + batch
     String LocalSampleDir = LocalAlignDir + "/" + SubDir
     String log = StagingDir + "log/" + Name + "_align.log"

     String outdir = OutputDir + "/" + SubDir
     String dragen_outdir = outdir + "/dragen"

     command {
         if [ ! -d "${LocalAlignDir}" ]; then
             /bin/mkdir ${LocalAlignDir}
         fi

         /bin/mkdir ${LocalSampleDir} && \
         /bin/mkdir ${outdir} && \
         /opt/edico/bin/dragen -r ${DragenRef} --tumor-fastq1 ${fastq1} --tumor-fastq2 ${fastq2} --RGSM-tumor ${SM} --RGID-tumor ${RG} --RGLB-tumor ${LB} \
         --umi-enable true --umi-min-supporting-reads ${readfamilysize}  --umi-enable-probability-model-merging=false \
         --umi-correction-scheme=random --umi-fuzzy-window-size=0 --umi-metrics-interval-file ${CoverageBed} \
         --qc-coverage-region-1 ${CoverageBed} --qc-coverage-reports-1 full_res --qc-coverage-ignore-overlaps true \
         --enable-map-align true --enable-sort true --enable-map-align-output true --gc-metrics-enable=true \
         --enable-variant-caller=true --vc-enable-umi-liquid true --vc-target-bed ${CoverageBed} --vc-somatic-hotspots ${Hotspot} \
         --vc-enable-triallelic-filter false --vc-combine-phased-variants-distance 3 \
         --enable-sv true --sv-call-regions-bed ${CoverageBed} --sv-exome true --sv-output-contigs true \
         --output-dir ${LocalSampleDir} --output-file-prefix ${Name} --output-format BAM &> ${log} && \
         /bin/mv ${log} ./ && \
         /bin/mv ${LocalSampleDir} ${dragen_outdir}
     }

     runtime {
         docker_image: DragenDockerImage
         dragen_env: DragenEnv
         cpu: "20"
         memory: "200 G"
         queue: queue
         job_group: jobGroup 
     }

     output {
         File bam = "${dragen_outdir}/${Name}_tumor.bam"
         File bai = "${dragen_outdir}/${Name}_tumor.bam.bai"
         File vcf = "${dragen_outdir}/${Name}.hard-filtered.vcf.gz"
         File index = "${dragen_outdir}/${Name}.hard-filtered.vcf.gz.tbi"
     }
}

task move_demux_fastq {
     input{
         Array[String] order_by
         String Batch
         String DemuxFastqDir
         String queue
         String jobGroup
     }

     String LocalDemuxFastqDir = "/staging/runs/MyeloSeqHD/demux_fastq/" + Batch

     command {
         if [ -d "${LocalDemuxFastqDir}" ]; then
             /bin/mv ${LocalDemuxFastqDir} ${DemuxFastqDir}
         fi
     }
     runtime {
         docker_image: "docker1(ubuntu:xenial)"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}

task batch_qc {
     input {
         Array[String] order_by
         String BatchDir
         String QC_pl
         String queue
         String jobGroup
     }
     command {
         if [ -n "$(/bin/ls -d ${BatchDir}/TW*)" ]; then
             /bin/chmod -R 777 ${BatchDir}/TW*
         fi

         /usr/bin/perl ${QC_pl} ${BatchDir}
     }
     runtime {
         docker_image: "docker1(mgibio/myeloseqhd:v2)"
         memory: "4 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}

task data_transfer {
     input {
         String order_by
         String BatchDir
         String xfer_pl
         String queue
         String jobGroup
         String? excluded
     }
     command {
         if [ -n "${excluded}" ]; then
             /usr/bin/perl ${xfer_pl} ${BatchDir} ${excluded}
         else
             /usr/bin/perl ${xfer_pl} ${BatchDir}
         fi
     }
     runtime {
         docker_image: "docker1(mgibio/data-transfer-helper:v1)"
         memory: "8 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}
