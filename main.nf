#!/usr/bin/env nextflow

/* 
 * Set up variables
 */

my_bed_ch = Channel
            .fromPath(params.input_bed)
            .ifEmpty { exit 1, "Cannot find input file : ${params.input_bed}" }

aggv2_bed_ch = Channel
            .fromPath(params.aggv2_chunks_bed)
            .ifEmpty { exit 1, "Cannot find input file : ${params.aggv2_chunks_bed}" }

Channel
      .fromPath(params.severity_scale)
      .ifEmpty { exit 1, "Cannot find severity scale : ${params.severity_scale}" }
      .set {severity_scale_ch}


/*
 * Start pipeline
 */

process find_chunk {
    
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file(my_bed) from my_bed_ch
    file(aggv2_bed) from aggv2_bed_ch

    output:
    file(geno_files) into geno_vcf_list_ch
    file(anno_files) into vep_vcf_list_ch

    shell:

    '''
    
	while read -r line; do
    gene="$(echo "${line}"| awk '{print $4}')";
	printf "${line}" > ${gene}.bed;
    gvcf="$(bedtools intersect -wo -a ${gene}.bed -b !{aggv2_bed} |cut -f 10)";
	gvcf_index="$(echo $(bedtools intersect -wo -a ${gene}.bed -b !{aggv2_bed} |cut -f 10).csi)";
    avcf="$(bedtools intersect -wo -a ${gene}.bed -b !{aggv2_bed} |cut -f 11)";
	avcf_index="$(echo $(bedtools intersect -wo -a ${gene}.bed -b !{aggv2_bed} |cut -f 11).csi)";
	echo "$gene,$gvcf,$gvcf_index" >> geno_files
	echo "$gene,$avcf,$avcf_index" >> anno_files
    done < !{my_bed}

    '''

}
/*
 * modify channel 
 */
vep_vcf_list_ch
		.splitCsv()
		.map {row -> tuple(val(row[0]), file(row[1]), file(row[2])) }
		.set {vep_vcf_ch}

process extract_variant_vep {

	publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(gene), file(avcf), file(avcf_index) from vep_vcf_ch
    file(severity_scale) from severity_scale_ch

    output:
    tuple file("${gene}_annotation.vcf.gz"), file("${gene}_annotation.vcf.gz.csi") into annotation_vcf_ch

    script:

    """
    bcftools +split-vep -i 'SYMBOL="'"${gene}"'"' -c SYMBOL -s worst:missense+ -S ${severity_scale} ${avcf} -O z -o ${gene}_annotation.vcf.gz
    bcftools index ${gene}_annotation.vcf.gz
    """

}
