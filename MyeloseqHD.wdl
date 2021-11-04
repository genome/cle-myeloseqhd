import "MyeloseqHDAnalysis.wdl" as subWF

workflow MyeloseqHD {

    File SampleSheet
    # sample sheet has this structure:
    # index  name  RG_ID  RG_FLOWCELL  RG_LANE  RG_LIB  RG_SAMPLE
    
    File DemuxSampleSheet

    Array[Array[String]] inputData = read_tsv(SampleSheet)
    Array[String] Adapters = ["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC","AGATCGGAAGAGCGTCGTGTAGGGAAA"]
    
    String IlluminaDir
    String JobGroup
    String OutputDir
    
    String Queue
    String DragenQueue = "duncavagee"
    
    String DragenReference = "/storage1/fs1/gtac-mgi/Active/CLE/reference/dragen_hg38"
    String Reference    = "/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/refdata/hg38/all_sequences.fa"
    String ReferenceDict = "/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/refdata/hg38/all_sequences.dict"
    
    String VEP          = "/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/VEP_cache"
    String QcMetrics    = "/storage1/fs1/gtac-mgi/Active/CLE/analysis/new_myeloseq/git/cle-myeloseqhd/accessory_files/MyeloseqHDQCMetrics.txt"
    String Description  = "/storage1/fs1/gtac-mgi/Active/CLE/analysis/new_myeloseq/git/cle-myeloseqhd/accessory_files/MyeloseqDescription.txt"
    
    String HaplotectBed = "/storage1/fs1/gtac-mgi/Active/CLE/analysis/new_myeloseq/git/cle-myeloseqhd/accessory_files/myeloseq.haplotect_snppairs_hg38.041718.bed"
    String AmpliconBed  = "/storage1/fs1/gtac-mgi/Active/CLE/analysis/new_myeloseq/git/cle-myeloseqhd/MyeloseqHD.16462-1615924889.Amplicons.hg38.110221.bed"
    String CoverageBed  = "/storage1/fs1/gtac-mgi/Active/CLE/analysis/new_myeloseq/git/cle-myeloseqhd/accessory_files/MyeloseqHD.16462-1615924889.CoverageQC.hg38.110221.bed"

    String CustomAnnotationVcf   = "/storage1/fs1/gtac-mgi/Active/CLE/analysis/new_myeloseq/git/cle-myeloseqhd/accessory_files/myeloseq_custom_annotations.annotated.011618.hg38.vcf.gz"
    String CustomAnnotationIndex = "/storage1/fs1/gtac-mgi/Active/CLE/analysis/new_myeloseq/git/cle-myeloseqhd/accessory_files/myeloseq_custom_annotations.annotated.011618.hg38.vcf.gz.tbi"

    String QC_pl = "/storage1/fs1/gtac-mgi/Active/CLE/analysis/new_myeloseq/git/cle-myeloseqhd/QC_metrics.pl"

    String CustomAnnotationParameters = "MYELOSEQ,vcf,exact,0,TCGA_AC,MDS_AC,MYELOSEQBLACKLIST"
    

    call dragen_demux {
        input: Dir=IlluminaDir,
               OutputDir=OutputDir,
               SampleSheet=DemuxSampleSheet,
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

    # the inputdata should be: index name id flow_cell lane lib sample_name
    scatter (samples in inputData){
        call get_fastq {
            input: Name=samples[1],
                   SampleSheetFile=prepare_samples.sample_sheet,
                   queue=Queue,
                   jobGroup=JobGroup
        }

        call dragen_align {
            input: DragenRef=DragenReference,
                   fastq1=get_fastq.demux_fastq1_file,
                   fastq2=get_fastq.demux_fastq2_file,
		   Name=samples[1],
		   RG=samples[3] + '.' + samples[4] + '.' + samples[0],
		   SM=samples[6],
		   LB=samples[5] + '.' + samples[0],
                   AmpliconBed=AmpliconBed,
                   CoverageBed=CoverageBed,
                   OutputDir=OutputDir,
                   SubDir=samples[1] + '_' + samples[0],
                   queue=DragenQueue,
                   jobGroup=JobGroup
        }

	call convert_bam {
	    input: Bam=dragen_align.bam,
             	   BamIndex=dragen_align.bai,
                   AmpliconBed=AmpliconBed,
                   refFasta=Reference,
                   Name=samples[1],
                   OutputDir=OutputDir,
                   SubDir=samples[1] + '_' + samples[0],
                   queue=Queue,
                   jobGroup=JobGroup
	}

        call subWF.MyeloseqHDAnalysis {
            input: Cram=convert_bam.cram,
                   CramIndex=convert_bam.crai,
                   CoverageBed=CoverageBed,
                   refFasta=Reference,
                   ReferenceDict=ReferenceDict,
                   Name=samples[1],
                   Vepcache=VEP,
                   CustomAnnotationVcf=CustomAnnotationVcf,
                   CustomAnnotationIndex=CustomAnnotationIndex,
                   CustomAnnotationParameters=CustomAnnotationParameters,
                   HaplotectBed=HaplotectBed,
                   QcMetrics=QcMetrics,
                   Description=Description,
                   OutputDir=OutputDir,
                   SubDir=samples[1] + '_' + samples[0],
                   Queue=Queue,
                   JobGroup=JobGroup
        }
    } 
}
        
#        call haloplex_qc {
#            input: refFasta=Reference,
#                   Bam=convert_bam.cram,
#                   AnnotatedTsv=run_vep.filtered_tsv,
#                   TargetBed=TargetBed,
#                   AmpliconBed=AmpliconBed,
#                   CoverageBed=CoverageBed,
#                   QcMetrics=QcMetrics,
#                   Description=Description,
#                   Name=samples[1],
#                   DemuxFile=prepare_samples.sample_sheet,
#                   Haplotect=run_haplotect.out,
#                   HaplotectSites=run_haplotect.sites,
#                   queue=Queue,
#                   jobGroup=JobGroup
#        }

#        call make_reports {
#            input: Variants=run_vep.filtered_tsv,
#                   QC=haloplex_qc.qc_json,
#                   Description=Description,
#                   Name=samples[1],
#                   queue=Queue,
#                   jobGroup=JobGroup
#        }

#        call igv_session {
#            input: Name=samples[1],
#                   queue=Queue,
#                   jobGroup=JobGroup
#        }

#   }

#    call batch_qc {
#        input: order_by=gather_files.out,
#               BatchDir=OutputDir,
#               QC_pl=QC_pl,
#               queue=Queue,
#               jobGroup=JobGroup
#    }


task dragen_demux {
     String Dir
     String OutputDir
     String SampleSheet 
     String jobGroup
     String queue

     String StagingDir = "/staging/runs/Haloplex/"
     String log = StagingDir + "log/dragen_demux.log"
     String DemuxReportDir = OutputDir + "/dragen_demux_reports"

     command <<<
          /opt/edico/bin/dragen --bcl-conversion-only true --bcl-only-matched-reads true --strict-mode true --sample-sheet ${SampleSheet} --bcl-input-directory ${Dir} --output-directory ./demux_fastq &> ${log} && \
          /bin/mv ${log} ./ && \
          /bin/mv ./demux_fastq/* ./ && \
          /bin/cp -r ./Reports ${DemuxReportDir}
     >>>

     runtime {
         docker_image: "seqfu/centos7-dragen-3.9.3:latest"
         cpu: "20"
         memory: "200 G"
         queue: queue
#        host: "compute1-dragen-3"
         job_group: jobGroup 
     }
     output {
         Array[File] read1 = glob("*_R1_001.fastq.gz") 
         Array[File] read2 = glob("*_R2_001.fastq.gz")
     }
}

task prepare_samples {
     File SampleSheet
     Array[File] Fastq1
     Array[File] Fastq2
     String jobGroup
     String queue

     command <<<
             /bin/cat ${write_tsv(Fastq1)} > 1.tmp.txt
             /bin/cat ${write_tsv(Fastq2)} > 2.tmp.txt
             /usr/bin/perl -e 'open(R1,"1.tmp.txt"); @r1 = <R1>; \
                 chomp @r1; close R1;\
                 open(R2,"2.tmp.txt"); @r2 = <R2>; \
                 chomp @r2; close R2; \
                 open(SS,"${SampleSheet}");
                 while(<SS>){
                     chomp;
                     my @l = split("\t",$_);
                     my $s = $l[1].'_';
                     my $r1 = (grep /$s/, @r1)[0];
                     my $r2 = (grep /$s/, @r2)[0];
                     my $persamplereads1 = `gunzip -c $r1 | wc -l`;
                     chomp $persamplereads1;
                     my $persamplereads2 = `gunzip -c $r2 | wc -l`;
                     chomp $persamplereads2;
                     print join("\t",$l[0],$l[1],$r1,$r2,($persamplereads1 / 4) + ($persamplereads2 / 4)),"\n";
                 }
                 close SS;' > sample_sheet.txt
     >>>
     runtime {
         docker_image: "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
         cpu: "1"
         memory: "4 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File sample_sheet = "sample_sheet.txt"
         Array[Array[String]] sample_data = read_tsv("sample_sheet.txt")
     }
}

task get_fastq {
     File SampleSheetFile
     String Name
     String jobGroup
     String queue

     String lookup = Name + "_"

     command {
         /bin/cp $(/bin/grep ${lookup} ${SampleSheetFile} | cut -f 3) "${Name}.1.fastq.gz" && \
         /bin/cp $(/bin/grep ${lookup} ${SampleSheetFile} | cut -f 4) "${Name}.2.fastq.gz"
     }
     runtime {
         docker_image: "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
         cpu: "1"
         memory: "8 G"
         queue: queue
         job_group: jobGroup 
     }
     output {
         File demux_fastq1_file = "${Name}.1.fastq.gz"
         File demux_fastq2_file = "${Name}.2.fastq.gz"
     }
}

task dragen_align {
     String Name
     String DragenRef
     String fastq1
     String fastq2
     String RG
     String SM
     String LB
     String AmpliconBed
     String CoverageBed
     Int? TrimLen
     String OutputDir
     String SubDir
     String jobGroup
     String queue

     String outdir = OutputDir + "/" + SubDir
     String dragen_outdir = outdir + "/dragen"
     String LocalSampleDir = "./" + SubDir

     String StagingDir = "/staging/runs/Haloplex/"
     String log = StagingDir + "log/" + Name + "_dragen_align.log"

     command {
         /bin/mkdir ${LocalSampleDir} && \
         /bin/mkdir ${outdir} && \
         /opt/edico/bin/dragen -r ${DragenRef} --tumor-fastq1 ${fastq1} --tumor-fastq2 ${fastq2} --RGSM-tumor ${SM} --RGID-tumor ${RG} --RGLB-tumor ${LB} --enable-map-align true --enable-sort true --enable-map-align-output true --vc-enable-umi-liquid true --gc-metrics-enable=true --qc-coverage-region-1 ${CoverageBed} --qc-coverage-reports-1 full_res --umi-enable true --umi-library-type=random-simplex --umi-min-supporting-reads 1 --enable-variant-caller=true --vc-target-bed ${CoverageBed} --umi-metrics-interval-file ${CoverageBed} --read-trimmers=fixed-len --trim-r1-5prime=${default=1 TrimLen} --trim-r1-3prime=${default=1 TrimLen} --trim-r2-5prime=${default=1 TrimLen} --trim-r2-3prime=${default=1 TrimLen} --output-dir ${LocalSampleDir} --output-file-prefix ${Name} --output-format BAM &> ${log} && \
         /bin/mv ${log} ./ && \
         /bin/mv ${LocalSampleDir} ${dragen_outdir}
     }
     
     runtime {
         docker_image: "seqfu/centos7-dragen-3.9.3:latest"
         cpu: "20"
         memory: "200 G"
         queue: queue
#        host: "compute1-dragen-3"
         job_group: jobGroup 
     }
     
     output {
         File bam = "${dragen_outdir}/${Name}_tumor.bam"
         File bai = "${dragen_outdir}/${Name}_tumor.bam.bai"
     }
}

task convert_bam {
     String Bam
     String BamIndex
     String Name
     String refFasta
     String AmpliconBed
     String OutputDir
     String SubDir
     String jobGroup
     String queue

     String outdir = OutputDir + "/" + SubDir

     command <<<
     	     /usr/local/bin/tagbam -v ${Bam} ${AmpliconBed} /tmp/tagged.bam > ${Name}.ampinfo.txt && \
	     /usr/local/bin/samtools view -T ${refFasta} -C -o "${Name}.cram" /tmp/tagged.bam && \
	     /usr/local/bin/samtools index "${Name}.cram" &&
	     (cut -f 5 ${Name}.ampinfo.txt && cut -f 4 ${AmpliconBed}) | sort | uniq -c | awk '!/\./ { print $2,$1-1; }' > ${Name}.ampcounts.txt && \
             /bin/cp ${Name}.ampinfo.txt ${outdir} && \
             /bin/cp ${Name}.ampcounts.txt ${outdir} && \
             /bin/cp ${Name}.cram ${outdir} && \
             /bin/cp ${Name}.cram.crai ${outdir}
     >>>
     
     runtime {
         docker_image: "registry.gsc.wustl.edu/mgi-cle/myeloseqhd:v1"
         cpu: "1"
         memory: "24 G"
         queue: queue
         job_group: jobGroup 
     }
     
     output {
         File cram  = "${Name}.cram"
         File crai = "${Name}.cram.crai"
	 File info = "${Name}.ampinfo.txt"
	 File counts = "${Name}.ampcounts.txt"
     }
}

task batch_qc {
     Array[String] order_by
     String BatchDir
     String QC_pl
     String queue
     String jobGroup

     command {
         /usr/bin/perl ${QC_pl} ${BatchDir}
     }  
     runtime {
         docker_image: "registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:compute1-20"
         memory: "4 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }   
}

