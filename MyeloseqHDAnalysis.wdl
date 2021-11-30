workflow MyeloseqHDAnalysis {

    String Cram
    String CramIndex
    String Name
    
    String refFasta 
    String Vepcache
    String CoverageBed
    String HaplotectBed
    String ReferenceDict

    String CustomAnnotationVcf
    String CustomAnnotationIndex
    String CustomAnnotationParameters

    String QcMetrics
    String Description

    String SubDir
    String OutputDir

    String Queue
    String JobGroup 

    call run_freebayes {
        input: Cram=Cram,
               CramIndex=CramIndex,
               CoverageBed=CoverageBed,
               refFasta=refFasta,
               Name=Name,
               queue=Queue,
               jobGroup=JobGroup
    }

    call clean_variants as clean_freebayes {
        input: Vcf=run_freebayes.vcf,
               Name=Name,
               refFasta=refFasta,
               queue=Queue,
               jobGroup=JobGroup
    }

    call run_pindel_region as run_pindel_flt3itd {
        input: Cram=Cram,
               CramIndex=CramIndex,
               Reg='chr13:28033987-28034316',
               refFasta=refFasta,
               Name=Name,
               queue=Queue,
               jobGroup=JobGroup
    }

    call bgzip_tabix as bgzip_tabix_pindel {
        input: Vcf=run_pindel_flt3itd.vcf,
               Name=Name,
               queue=Queue,
               jobGroup=JobGroup
    }

    call clean_variants as clean_pindel_itd {
        input: Vcf=bgzip_tabix_pindel.vcf,
               Name=Name,
               refFasta=refFasta,
               queue=Queue,
               jobGroup=JobGroup
    }

    call combine_variants {
        input: Vcfs=[clean_freebayes.cleaned_vcf_file,clean_pindel_itd.cleaned_vcf_file],
               Cram=Cram,
               CramIndex=CramIndex,
               refFasta=refFasta,
               Name=Name,
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

    call run_haplotect {
        input: refFasta=refFasta,
               refDict=ReferenceDict,
               Cram=Cram,
               CramIndex=CramIndex,
               Bed=HaplotectBed,
               Name=Name,
               queue=Queue,
               jobGroup=JobGroup
    }

    call gather_files {
        input: OutputFiles=[clean_pindel_itd.cleaned_vcf_file,
               combine_variants.vcf,
               run_haplotect.out_file,
               run_haplotect.sites_file,
               run_vep.vcf,
               run_vep.filtered_vcf,
               run_vep.filtered_tsv],
               OutputDir=OutputDir,
               SubDir=SubDir,
               queue=Queue,
               jobGroup=JobGroup
    }

    call haloplex_qc {
      input: order_by=gather_files.done,
      refFasta=refFasta,
      Name=Name,
      CoverageBed=CoverageBed,
      QcMetrics=QcMetrics,
      Description=Description,
      OutputDir=OutputDir,
      SubDir=SubDir,
      queue=Queue,
      jobGroup=JobGroup
    }

    output {
        String all_done = haloplex_qc.done
    }
}

task run_freebayes {
     File Cram
     File CramIndex
     String CoverageBed
     String refFasta
     String Name
     Float? MinFreq
     Int? MinReads
     Int? MinMapQual
     Int? MaxMismatch

     String jobGroup
     String queue

     command {
       /usr/local/bin/freebayes -C ${default=3 MinReads} -q ${default=13 MinMapQual} -F ${default="0.0008" MinFreq} -$ ${default=4 MaxMismatch} \
       -f ${refFasta} -t ${CoverageBed} ${Cram} > "${Name}.freebayes.vcf"
     }

     runtime {
         docker_image: "registry.gsc.wustl.edu/mgi-cle/myeloseqhd:v1"
         cpu: "1"
         memory: "16 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File vcf = "${Name}.freebayes.vcf"
     }
}

task run_pindel_region {
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

    command <<<
        (set -eo pipefail && /usr/local/bin/samtools view -T ${refFasta} ${Cram} ${Reg} | /opt/pindel-0.2.5b8/sam2pindel - /tmp/in.pindel ${default=250 Isize} tumor 0 Illumina-PairEnd) && \
        /usr/bin/pindel -f ${refFasta} -p /tmp/in.pindel -c ${Reg} -o /tmp/out.pindel && \
        /usr/bin/pindel2vcf -P /tmp/out.pindel -G -r ${refFasta} -e ${default=3 MinReads} -R ${default="GRCh38" Genome} -d ${default="GRCh38" Genome} -v /tmp/out.vcf && \
        sed 's/END=[0-9]*;//' /tmp/out.vcf > ${Name}.pindel.vcf
    >>>

    runtime {
        docker_image: "registry.gsc.wustl.edu/fdu/pindel2vcf-0.6.3:1"
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
     String Vcf
     String Name
     String queue
     String jobGroup

     command {
          /opt/htslib/bin/bgzip -c ${Vcf} > ${Name}.bgzip_tabix.vcf.gz && \
          /usr/bin/tabix -p vcf ${Name}.bgzip_tabix.vcf.gz
     }
     runtime {
         docker_image: "registry.gsc.wustl.edu/mgi-cle/myeloseqhd:v1"
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
     String Vcf
     String Name
     String refFasta
     String jobGroup
     String queue

     command {
         /usr/local/bin/bcftools sort -Oz ${Vcf} | \
         /usr/local/bin/bcftools norm -m-any -f ${refFasta} -Oz | \
         /usr/local/bin/bcftools norm -d any -Oz > "${Name}.cleaned.vcf.gz" && \
         /usr/bin/tabix -p vcf "${Name}.cleaned.vcf.gz"
     }

     runtime {
         docker_image: "registry.gsc.wustl.edu/mgi-cle/myeloseqhd:v1"
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
     Array[String] Vcfs
     String Cram
     String CramIndex
     String refFasta
     String Name
     String jobGroup
     String queue

     command {
         /usr/local/bin/bcftools merge --force-samples -Oz -o combined.vcf.gz ${sep=" " Vcfs} && \
         /usr/bin/tabix -p vcf combined.vcf.gz && \
         /usr/bin/python3 /usr/local/bin/filterHaloplex.py -r ${refFasta} combined.vcf.gz ${Cram} ${Name} > ${Name}.combined_and_tagged.vcf
     }

     runtime {
         docker_image: "registry.gsc.wustl.edu/mgi-cle/myeloseqhd:v1"
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

     command {
         if [ $(/bin/grep -v '^#' ${CombineVcf}|/usr/bin/wc -l) == 0 ]; then
             /bin/cp ${CombineVcf} ${Name}.annotated.vcf && \
             /bin/cp ${CombineVcf} ${Name}.annotated_filtered.vcf && \
             /opt/htslib/bin/bgzip ${Name}.annotated.vcf && /usr/bin/tabix -p vcf ${Name}.annotated.vcf.gz && \
             /opt/htslib/bin/bgzip ${Name}.annotated_filtered.vcf && /usr/bin/tabix -p vcf ${Name}.annotated_filtered.vcf.gz && \
             /usr/bin/touch ${Name}.variants_annotated.tsv
         else
             /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl --format vcf \
             --vcf --plugin Downstream --plugin Wildtype --fasta ${refFasta} --hgvs --symbol --term SO --flag_pick \
             -i ${CombineVcf} --custom ${CustomAnnotationVcf},${CustomAnnotationParameters} --offline --cache --max_af --dir ${Vepcache} -o ${Name}.annotated.vcf && \
             /opt/htslib/bin/bgzip ${Name}.annotated.vcf && /usr/bin/tabix -p vcf ${Name}.annotated.vcf.gz && \
             /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /opt/vep/ensembl-vep/filter_vep -i ${Name}.annotated.vcf.gz --format vcf \
             --filter "(MAX_AF < ${default='0.001' maxAF} or not MAX_AF) or MYELOSEQ_TCGA_AC or MYELOSEQ_MDS_AC" -o ${Name}.annotated_filtered.vcf && \
             /opt/htslib/bin/bgzip ${Name}.annotated_filtered.vcf && /usr/bin/tabix -p vcf ${Name}.annotated_filtered.vcf.gz && \
             /usr/bin/java -Xmx4g -jar /opt/GenomeAnalysisTK.jar -T VariantsToTable \
             -R ${refFasta} --showFiltered --variant ${Name}.annotated_filtered.vcf.gz -o /tmp/variants.tsv \
             -F CHROM -F POS -F ID -F FILTER -F REF -F ALT -GF TAMP -GF SAMP -GF CVAF -GF RO -GF AO && \
             /usr/bin/python /usr/bin/add_annotations_to_table_helper.py /tmp/variants.tsv ${Name}.annotated_filtered.vcf.gz \
             Consequence,SYMBOL,EXON,INTRON,Feature_type,Feature,HGVSc,HGVSp,HGNC_ID,MAX_AF,MYELOSEQ_TCGA_AC,MYELOSEQ_MDS_AC /tmp/ && \
             mv /tmp/variants.annotated.tsv ${Name}.variants_annotated.tsv
         fi
     }
     runtime {
         docker_image: "registry.gsc.wustl.edu/fdu/vep90-gatk3.6-htslib1.3.2:1"
         cpu: "1"
         memory: "10 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File vcf = "${Name}.annotated.vcf.gz"
         File filtered_vcf = "${Name}.annotated_filtered.vcf.gz"
         File filtered_tsv = "${Name}.variants_annotated.tsv"
     }
}

task run_haplotect {
     String Cram
     String CramIndex
     String Bed
     String Name
     String refDict
     String refFasta
     String queue
     String jobGroup

     Int? MinReads

     command <<<
             /usr/bin/awk -v OFS="\t" '{ $2=$2-1; print; }' ${Bed} > /tmp/pos.bed && \
             /usr/local/openjdk-8/bin/java -Xmx6g \
             -jar /opt/hall-lab/gatk-package-4.1.8.1-18-ge2f02f1-SNAPSHOT-local.jar Haplotect \
             -I ${Cram} -R ${refFasta} --sequence-dictionary ${refDict} \
             -mmq 20 -mbq 20 -max-depth-per-sample 10000 -gstol 0.001 -mr ${default=10 MinReads} \
             -htp ${Bed} -L /tmp/pos.bed -outPrefix ${Name}
     >>>

     runtime {
             docker_image: "registry.gsc.wustl.edu/mgi-cle/haplotect:0.3"
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

task haloplex_qc {
  String order_by
  String refFasta
  String Name
  String CoverageBed
  String QcMetrics
  String Description
  String OutputDir
  String SubDir
  String jobGroup
  String queue

  String SampleOutDir = OutputDir + "/" + SubDir

  command {
    /usr/bin/perl /usr/local/bin/CalculateCoverageQC.pl -r ${refFasta} -d ${SampleOutDir} -n ${Name} \
    -t ${CoverageBed} -q ${QcMetrics} -i ${Description} && \
    /bin/mv ./*.qc.txt ./*.qc.json ${SampleOutDir}
  }
  runtime {
    docker_image: "registry.gsc.wustl.edu/mgi-cle/myeloseqhd:v1"
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
     Array[String] OutputFiles
     String OutputDir
     String? SubDir
     String jobGroup
     String queue

     command {
         if [[ ${SubDir} != "" ]] && [[ ! -e ${OutputDir}/${SubDir} ]]; then
             mkdir ${OutputDir}/${SubDir}
         fi
         /bin/mv -f -t ${OutputDir}/${SubDir} ${sep=" " OutputFiles}
     }
     runtime {
         docker_image: "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}
