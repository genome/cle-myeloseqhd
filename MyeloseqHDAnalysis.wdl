version 1.0

workflow MyeloseqHDAnalysis {

    input {

        String Bam
        String BamIndex
        String DragenVcf
        String DragenVcfIndex
        String Name

        String refFasta 
        String ReferenceDict
        String Vepcache
        String VariantDB

        String MyeloSeqHDRepo

        String HaplotectBed = MyeloSeqHDRepo + "/accessory_files/MyeloseqHD.haplotect.bed"
        String AmpliconBed  = MyeloSeqHDRepo + "/accessory_files/MyeloseqHD.amplicons.bed"
        String CoverageBed  = MyeloSeqHDRepo + "/accessory_files/MyeloseqHD.CoverageQC.hg38.bed"
        String GenotypeVcf  = MyeloSeqHDRepo + "/accessory_files/MyeloseqHD.forcegenotype.vcf.gz"
        String QcMetrics    = MyeloSeqHDRepo + "/accessory_files/MyeloseqHD.QCMetrics.json"
        String CustomAnnotationVcf   = MyeloSeqHDRepo + "/accessory_files/MyeloseqHD.custom_annotations.vcf.gz"
        String CustomAnnotationIndex = MyeloSeqHDRepo + "/accessory_files/MyeloseqHD.custom_annotations.vcf.gz.tbi"
        String CustomAnnotationParameters = "MYELOSEQ,vcf,exact,0,TCGA_AC,MDS_AC,BLACKLIST"

        Int MinReads
        Float MinVaf

        String mrn
        String all_mrn
        String accession
        String DOB
        String sex
        String exception
        String RunInfoString

        String SubDir
        String OutputDir

        String Queue
        String JobGroup 

    }

    call query_DB {
        input: mrn=all_mrn,
               accession=accession,
               VariantDB=VariantDB,
               queue=Queue,
               jobGroup=JobGroup
    }

    call clean_variants as clean_queryDB_vcf {
        input: Vcf=query_DB.query_vcf,
               Name=Name,
               refFasta=refFasta,
               queue=Queue,
               jobGroup=JobGroup
    }

    call convert_bam {
        input: Bam=Bam,
               BamIndex=BamIndex,
               AmpliconBed=AmpliconBed,
               refFasta=refFasta,
               Name=Name,
               OutputDir=OutputDir,
               SubDir=SubDir,
               queue=Queue,
               jobGroup=JobGroup
    }

    call run_pindel_region {
        input: Cram=convert_bam.cram,
               CramIndex=convert_bam.crai,
               Reg=CoverageBed,
               refFasta=refFasta,
               Name=Name,
               queue=Queue,
               jobGroup=JobGroup
    }

    call bgzip_tabix as bgzip_tabix_pindel {
        input: Vcf=run_pindel_region.vcf,
               Name=Name,
               queue=Queue,
               jobGroup=JobGroup
    }

    call clean_variants as clean_pindel_region {
        input: Vcf=bgzip_tabix_pindel.vcf,
               Name=Name,
               refFasta=refFasta,
               queue=Queue,
               jobGroup=JobGroup
    }

    call combine_variants {
        input: Vcfs=[DragenVcf,clean_pindel_region.cleaned_vcf_file,clean_queryDB_vcf.cleaned_vcf_file,GenotypeVcf],
               Cram=convert_bam.cram,
               CramIndex=convert_bam.crai,
               refFasta=refFasta,
               Vaf=MinVaf,
               Reads=MinReads,
               Name=Name,
               MyeloSeqHDRepo=MyeloSeqHDRepo,
               queue=Queue,
               jobGroup=JobGroup
    }

    call run_vep {
        input: CombineVcf=combine_variants.vcf,
               refFasta=refFasta,
               Vepcache=Vepcache,
               CustomAnnotationVcf=CustomAnnotationVcf,
               CustomAnnotationIndex=CustomAnnotationIndex,
               CustomAnnotationParameters=CustomAnnotationParameters,
               Name=Name,
               queue=Queue,
               jobGroup=JobGroup
    }

    call upload_DB {
        input: Vcf=run_vep.filtered_vcf,
               mrn=mrn,
               accession=accession,
               VariantDB=VariantDB,
               CoverageBed=CoverageBed,
               queue=Queue,
               jobGroup=JobGroup
    }

    call run_haplotect {
        input: refFasta=refFasta,
               refDict=ReferenceDict,
               Cram=convert_bam.cram,
               CramIndex=convert_bam.crai,
               Bed=HaplotectBed,
               Name=Name,
               queue=Queue,
               jobGroup=JobGroup
    }

    call gather_files {
        input: order_by=upload_DB.done,
               OutputFiles=[run_haplotect.out_file,
               run_haplotect.sites_file,
               run_vep.vcf,
               run_vep.filtered_vcf],
               OutputDir=OutputDir,
               SubDir=SubDir,
               queue=Queue,
               jobGroup=JobGroup
    }

    call make_report {
        input: order_by=gather_files.done,
               Name=Name,
               mrn=mrn,
               accession=accession,
               DOB=DOB,
               sex=sex,
               exception=exception,
               RunInfoString=RunInfoString,
               CoverageBed=CoverageBed,
               QcMetrics=QcMetrics,
               MyeloSeqHDRepo=MyeloSeqHDRepo,
               OutputDir=OutputDir,
               SubDir=SubDir,
               queue=Queue,
               jobGroup=JobGroup
    }

    output {
        String all_done = make_report.done
    }
}


task convert_bam {

    input {
        String Bam
        String BamIndex
        String Name
        String refFasta
        String AmpliconBed
        String OutputDir
        String SubDir
        String jobGroup
        String queue
    }

    String outdir = OutputDir + "/" + SubDir

     command <<<
         /usr/local/bin/tagbam -v ~{Bam} ~{AmpliconBed} /tmp/tagged.bam > ~{Name}.ampinfo.txt && \
         /usr/local/bin/samtools view -T ~{refFasta} -C -o "~{Name}.cram" /tmp/tagged.bam && \
         /usr/local/bin/samtools index "~{Name}.cram" &&
         (cut -f 5 ~{Name}.ampinfo.txt && cut -f 4 ~{AmpliconBed}) | sort | uniq -c | awk '!/\./ { print $2,$1-1; }' > ~{Name}.ampcounts.txt && \
         /bin/cp ~{Name}.ampinfo.txt ~{outdir} && \
         /bin/cp ~{Name}.ampcounts.txt ~{outdir} && \
         /bin/cp ~{Name}.cram ~{outdir} && \
         /bin/cp ~{Name}.cram.crai ~{outdir}
     >>>

     runtime {
         docker_image: "docker1(mgibio/myeloseqhd:v2)"
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

task run_pindel_region {

    input {
        File Cram
        File CramIndex
        String Reg
        String? Genome
        Int? Isize
        Int? MinReads
        String refFasta
        String Name
        String jobGroup
        String queue
    }

     command <<<
         (set -eo pipefail && /usr/local/bin/samtools view -T ~{refFasta} -ML ~{Reg} ~{Cram} | /opt/pindel-0.2.5b8/sam2pindel - /tmp/in.pindel ~{default=250 Isize} tumor 0 Illumina-PairEnd) && \
         /usr/local/bin/pindel -f ~{refFasta} -p /tmp/in.pindel -j ~{Reg} -o /tmp/out.pindel && \
         /usr/local/bin/pindel2vcf -P /tmp/out.pindel -G -r ~{refFasta} -e ~{default=3 MinReads} -R ~{default="GRCh38" Genome} -d ~{default="GRCh38" Genome} -v /tmp/out.vcf && \
         sed 's/END=[0-9]*;//' /tmp/out.vcf > ~{Name}.pindel.vcf
     >>>

     runtime {
         docker_image: "docker1(mgibio/pindel2vcf-0.6.3:1)"
         cpu: "1"
         memory: "16 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File vcf = "${Name}.pindel.vcf"
     }
}

task bgzip_tabix {
    input{
        String Vcf
        String Name
        String queue
        String jobGroup
    }

     command {
         /opt/htslib/bin/bgzip -c ${Vcf} > ${Name}.bgzip_tabix.vcf.gz && \
         /usr/bin/tabix -p vcf ${Name}.bgzip_tabix.vcf.gz
     }
     runtime {
         docker_image: "docker1(mgibio/myeloseqhd:v2)"
         cpu: "1"
         memory: "8 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File vcf = "${Name}.bgzip_tabix.vcf.gz"
         File index = "${Name}.bgzip_tabix.vcf.gz.tbi"
     }
}

task clean_variants {
    input{
        String Vcf
        String Name
        String refFasta
        String jobGroup
        String queue
    }

     command {
         /usr/local/bin/bcftools sort -Oz ${Vcf} | \
         /usr/local/bin/bcftools norm -m-any -f ${refFasta} -Oz | \
         /usr/local/bin/bcftools norm -d any -Oz > "${Name}.cleaned.vcf.gz" && \
         /usr/bin/tabix -p vcf "${Name}.cleaned.vcf.gz"
     }

     runtime {
         docker_image: "docker1(mgibio/myeloseqhd:v2)"
         cpu: "1"
         memory: "16 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File cleaned_vcf_file = "${Name}.cleaned.vcf.gz"
         File index = "${Name}.cleaned.vcf.gz.tbi"
     }
}

task combine_variants {
    input{
        Array[String] Vcfs
        String MyeloSeqHDRepo
        String Cram
        String CramIndex
        String refFasta
        String Name
        String jobGroup
        String queue

        Int? Reads
        Int? MinReadsPerFamily
        Float? Vaf
    }

     command {
         /usr/local/bin/bcftools merge -F x -m none --force-samples -Oz ${sep=" " Vcfs} | /usr/local/bin/bcftools sort -Oz -o combined.vcf.gz && \
         /usr/bin/tabix -p vcf combined.vcf.gz && \
         /usr/bin/python3 ${MyeloSeqHDRepo}/dockerfiles/docker-myeloseq/filterHaloplex.py -r ${refFasta} --minreadsperfamily ${default='3' MinReadsPerFamily} -m ${default='5' Reads} -d ${default='0.02' Vaf} combined.vcf.gz ${Cram} ${Name} > ${Name}.combined_and_tagged.vcf
     }

     runtime {
         docker_image: "docker1(mgibio/myeloseqhd:v2)"
         cpu: "1"
         memory: "10 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File vcf = "${Name}.combined_and_tagged.vcf"
     }
}

task run_vep {
    input{
        File CombineVcf
        String refFasta
        String Vepcache
        File CustomAnnotationVcf
        File CustomAnnotationIndex
        String CustomAnnotationParameters
        Float? maxAF
        String Name
        String jobGroup
        String queue
    }

     command {
         if [ $(/bin/grep -v '^#' ${CombineVcf}|/usr/bin/wc -l) == 0 ]; then
             /bin/cp ${CombineVcf} ${Name}.annotated.vcf && \
             /bin/cp ${CombineVcf} ${Name}.annotated_filtered.vcf && \
             /usr/local/bin/bgzip ${Name}.annotated.vcf && /usr/local/bin/tabix -p vcf ${Name}.annotated.vcf.gz && \
             /usr/local/bin/bgzip ${Name}.annotated_filtered.vcf && /usr/local/bin/tabix -p vcf ${Name}.annotated_filtered.vcf.gz && \
             /usr/bin/touch ${Name}.variants_annotated.tsv
         else
             /usr/bin/perl -I /opt/vep/lib/perl/VEP/Plugins /opt/vep/src/ensembl-vep/vep --format vcf \
             --vcf --plugin Downstream --fasta ${refFasta} --hgvs --symbol --term SO --flag_pick \
             -i ${CombineVcf} --custom ${CustomAnnotationVcf},${CustomAnnotationParameters} --offline --cache --max_af --dir ${Vepcache} -o ${Name}.annotated.vcf && \
             /usr/local/bin/bgzip ${Name}.annotated.vcf && /usr/local/bin/tabix -p vcf ${Name}.annotated.vcf.gz && \
             /usr/bin/perl -I /opt/vep/lib/perl/VEP/Plugins /opt/vep/src/ensembl-vep/filter_vep -i ${Name}.annotated.vcf.gz --format vcf \
             --filter "(MAX_AF < ${default='0.001' maxAF} or not MAX_AF) or MYELOSEQ_TCGA_AC or MYELOSEQ_MDS_AC" -o ${Name}.annotated_filtered.vcf && \
             /usr/local/bin/bgzip ${Name}.annotated_filtered.vcf && /usr/local/bin/tabix -p vcf ${Name}.annotated_filtered.vcf.gz
         fi
     }
     runtime {
         docker_image: "docker1(mgibio/vep105-htslib1.9:1)"
         cpu: "1"
         memory: "10 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File vcf = "${Name}.annotated.vcf.gz"
         File filtered_vcf = "${Name}.annotated_filtered.vcf.gz"
     }
}

task run_haplotect {
    input{
        String Cram
        String CramIndex
        String Bed
        String Name
        String refDict
        String refFasta
        String queue
        String jobGroup

        Int? MinReads
    }

    command <<<
         /usr/bin/awk -v OFS="\t" '{ $2=$2-1; print; }' ~{Bed} > /tmp/pos.bed && \
         /usr/local/openjdk-8/bin/java -Xmx6g \
         -jar /opt/hall-lab/gatk-package-4.1.8.1-18-ge2f02f1-SNAPSHOT-local.jar Haplotect \
         -I ~{Cram} -R ~{refFasta} --sequence-dictionary ~{refDict} \
         -mmq 20 -mbq 20 -max-depth-per-sample 10000 -gstol 0.001 -mr ~{default=10 MinReads} \
         -htp ~{Bed} -L /tmp/pos.bed -outPrefix ~{Name}
    >>>

     runtime {
         docker_image: "docker1(abelhj/haplotect:0.3)"
         cpu: "1"
         memory: "8 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File out_file = "${Name}.haplotect.txt"
         File sites_file = "${Name}.haplotectloci.txt"
     }
}

task make_report {
    input{
        String order_by
        String Name
        String mrn
        String accession
        String DOB
        String sex
        String exception
        String RunInfoString
        String CoverageBed
        String QcMetrics
        String OutputDir
        String SubDir
        String MyeloSeqHDRepo
        String jobGroup
        String queue

        Int? MinReadsPerFamily
    }

     String SampleOutDir = OutputDir + "/" + SubDir

     command {
         /usr/bin/python3 ${MyeloSeqHDRepo}/dockerfiles/docker-myeloseq/make_hd_report.py -n ${Name} -d ${SampleOutDir} -c ${CoverageBed} -q ${QcMetrics} \
         -m ${mrn} -a ${accession} -b ${DOB} -e ${exception} -i ${RunInfoString} && \
         /bin/mv ./*.report.txt ./*.report.json ${SampleOutDir}
     }
     runtime {
         docker_image: "docker1(mgibio/myeloseqhd:v2)"
         cpu: "1"
         memory: "16 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}

task gather_files {
    input{
        String order_by
        Array[String] OutputFiles
        String OutputDir
        String? SubDir
        String jobGroup
        String queue
    }

     command {
         if [[ ${SubDir} != "" ]] && [[ ! -e ${OutputDir}/${SubDir} ]]; then
             mkdir ${OutputDir}/${SubDir}
         fi
         /bin/mv -f -t ${OutputDir}/${SubDir} ${sep=" " OutputFiles}
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

task upload_DB {
    input{
        String Vcf
        String mrn
        String accession
        String VariantDB
        String CoverageBed
        String jobGroup
        String queue
    }

     command {
         /usr/bin/python3 /usr/local/bin/variantDB.py -d ${VariantDB} -v ${Vcf} -c ${CoverageBed} -i ${mrn} -j ${accession}
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

task query_DB {
    input{
        String mrn
        String accession
        String VariantDB
        String jobGroup
        String queue
    }

     command {
         /usr/bin/python3 /usr/local/bin/variantDB.py -d ${VariantDB} -m ${mrn} -a ${accession}
     }
     runtime {
         docker_image: "docker1(mgibio/myeloseqhd:v2)"
         memory: "4 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File query_vcf = "${accession}_query.vcf"
     }
}
