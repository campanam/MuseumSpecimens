picard = "picard " + params.picard_java
gatk = 'gatk --java-options "' + params.gatk_java + '" '

process prepareRef {

	// Prepare reference sequence for downstream processing
		
	input:
	path refseq
	
	output:
	path "${refseq.baseName}*.{amb,ann,bwt,pac,sa,fai,dict}"
	
	"""
	bwa index ${refseq}
	samtools faidx ${refseq}
	samtools dict ${refseq} > ${refseq.baseName}.dict
	"""

}

process genMapIndex {

	// Generate GenMap index
	
	label 'genmap'
		
	input:
	path refseq
	val gm_tmpdir
	
	output:
	tuple path("$refseq"), path("${refseq.simpleName}_index"), path("${refseq.simpleName}_index/*")
	
	"""
	export TMPDIR=${gm_tmpdir}
	if [ ! -d ${gm_tmpdir} ]; then mkdir ${gm_tmpdir}; fi
	genmap index -F ${refseq} -I ${refseq.simpleName}_index
	"""

}

process genMapMap {

	// Calculate mappability using GenMap and filter using filterGM
	
	label 'genmap'
	label 'ruby'
		
	input:
	tuple path(refseq), path(genmap_index), path("*")
	
	output:
	path "${refseq.simpleName}_genmap.1.0.bed"
	
	"""
	genmap map ${params.gm_opts} -T ${task.cpus} -I ${refseq.simpleName}_index/ -O ${refseq.simpleName}_genmap -b
	filterGM.rb ${refseq.simpleName}_genmap.bed 1.0 exclude > ${refseq.simpleName}_genmap.1.0.bed
	"""
}

process trimPEAdapters {

	// Trim adapters for paired-end libraries using AdapterRemoval 2.3.3

	input:
	tuple val(sample), val(library), path(reads1), path(reads2), val(adapter1), val(adapter2), val(rg)
	
	output:
	tuple val(sample), val(library), path("${library}.fastq.gz"), val(rg)
	
	"""
	AdapterRemoval --file1 $reads1 --file2 $reads2 --basename $library --adapter1 $adapter1 --adapter2 $adapter2 --gzip --collapse --minlength 30
	cat ${library}.collapsed.gz ${library}.collapsed.truncated.gz > ${library}.fastq.gz
	"""

}

process trimSEAdapters {

	// Trim adapters for single-end libraries using AdapterRemoval 2.3.3

	input:
	tuple val(sample), val(library), path(reads1), val(adapter1), val(adapter2), val(rg)
	
	output:
	tuple val(sample), val(library), path("${library}.truncated.gz"), val(rg)
	
	"""
	AdapterRemoval --file1 $reads1 --basename $library --adapter1 $adapter1 --adapter2 $adapter2 --gzip --minlength 30
	"""

}

process alignSeqs {

	// Align fastqs against reference sequence
		
	input:
	tuple val(sample), val(library), path(reads1), val(rg)
	path refseq
	path "*"
	
	output:
	path("${library}_${refseq.simpleName}.bam"), emit: bam
	val(sample), emit: sample
	
	script:
	samtools_extra_threads = task.cpus - 1
	"""
	bwa samse -r '${rg}' ${refseq} <(bwa aln -t ${task.cpus} -l 1024 ${refseq} ${reads1}) ${reads1} | samtools fixmate -@ ${samtools_extra_threads} -m - - | samtools sort -@ ${samtools_extra_threads} -o ${library}_${refseq.simpleName}.bam - 
	"""
	
}

process realignIndels {

	// GATK RealignerTargetCreator. Index bam with picard.
		
	input:
	path rg_bam
	val sample
	path refseq
	path "*"
	
	output:
	tuple path("${rg_bam.simpleName}.realn.bam"), val(sample)
	
	"""
	$picard BuildBamIndex I=${rg_bam}
	$gatk LeftAlignIndels -R ${refseq} -I $rg_bam -O ${rg_bam.simpleName}.realn.bam
	"""

}

process markDuplicates {

	// Mark duplicates using sambamba, samtools or picard
	
	publishDir  "$params.outdir/01_OriginalBAMs", mode: 'copy'
	
	input:
	tuple path(sorted_bam), val(sample)
	
	output:
	tuple path("${sorted_bam.simpleName}.markdup.bam"), val(sample)
	
	script:
	samtools_extra_threads = task.cpus - 1
	if ( params.markDuplicates == "sambamba" )
		"""
		sambamba markdup ${sorted_bam} ${sorted_bam.simpleName}.markdup.bam
		"""
	else if ( params.markDuplicates == "samtools" )
		"""
		samtools markdup -@ ${samtools_extra_threads} ${sorted_bam} ${sorted_bam.simpleName}.markdup.bam
		"""
	else
		"""
		$picard MarkDuplicates I=${sorted_bam} O=${sorted_bam.simpleName}.markdup.bam M=${sorted_bam.simpleName}.markdup.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
		"""
		
}

process profileDamage {

	// Profile deamination damage for ancient libraries
	
	publishDir "$params.outdir/02_DamageProfiles", mode: 'copy', pattern: "*_damage/*.*"
	
	input:
	tuple path(mrkdupbam), val(sample)
	
	output:
	path "${mrkdupbam.simpleName}_damage/*pdf"
	path "${mrkdupbam.simpleName}_damage/*txt"
	path "${mrkdupbam.simpleName}_damage/*log"	
	"""
	damageprofiler ${params.java11_options} -i ${mrkdupbam} -o ${mrkdupbam.simpleName}_damage -r ${params.refseq}
	"""

}

process mergeLibraries {

	// Merge libraries by their sample IDs using SAMtools merge
	
	input:
	tuple path(bam), val(sample)
	
	output:
	path "${sample}_merged.bam"
	
	script:
	samtools_extra_threads = task.cpus - 1
	bams = 0
	bamlist = ""
	for (i in bam) {
		bams++
		bamlist = bamlist + " " + i
	}
	if (bams == 1) // Skip merging single libraries
		"""
		ln -s $bamlist ${sample}_merged.bam
		"""
	else
		"""
		samtools merge -@ ${samtools_extra_threads} -o ${sample}_merged.bam $bamlist
		"""
} 

process reRealignIndels {

	// GATK RealignerTargetCreator. Index bam with picard.
		
	input:
	path rg_bam
	path refseq
	path "*"
	
	output:
	path("${rg_bam.simpleName}.realn.bam")
	
	"""
	$picard BuildBamIndex I=${rg_bam}
	$gatk LeftAlignIndels -R ${refseq} -I $rg_bam -O ${rg_bam.simpleName}.realn.bam
	"""

}

process reMarkDuplicates {

	// Mark duplicates using sambamba, samtools or picard
	
	publishDir  "$params.outdir/03_MergedBAMs", mode: 'copy'
	
	input:
	path(sorted_bam)
	
	output:
	path("${sorted_bam.simpleName}.markdup.bam")
	
	script:
	samtools_extra_threads = task.cpus - 1
	if ( params.markDuplicates == "sambamba" )
		"""
		sambamba markdup ${sorted_bam} ${sorted_bam.simpleName}.markdup.bam
		"""
	else if ( params.markDuplicates == "samtools" )
		"""
		samtools markdup -@ ${samtools_extra_threads} ${sorted_bam} ${sorted_bam.simpleName}.markdup.bam
		"""
	else
		"""
		$picard MarkDuplicates I=${sorted_bam} O=${sorted_bam.simpleName}.markdup.bam M=${sorted_bam.simpleName}.markdup.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
		"""
		
}

process trimAncientTermini {

	// Trim a set number of bases from termini of ancient samples using bamUtil trimBam
	// Could theoretically calculate this from the DamageProfiler output, but a hard cut-off should work well enough
	
	publishDir "$params.outdir/04_AncientLibraryTrimmedBAMs", mode: 'copy'
	
	input:
	path(profbam)
	
	output:
	path("${profbam.simpleName}.trim.bam")
	
	"""
	bam trimBam $profbam ${profbam.simpleName}.trim.bam ${params.aDNA_trimmed_bases}
	"""
	
}

process calculateStatistics {

	// Calculate alignment statistics, depth, and coverage statistics using SAMtools
	
	publishDir "$params.outdir/05_AlignmentStatistics", mode: 'copy'
	
	input:
	path(trimbam)
	
	output:
	path("${trimbam.simpleName}.*.txt")
	
	script:
	samtools_extra_threads = task.cpus - 1
	"""
	samtools flagstat -@ ${samtools_extra_threads} $trimbam > ${trimbam.simpleName}.stats.txt
	samtools depth -@ ${samtools_extra_threads} $trimbam > ${trimbam.simpleName}.depth.txt
	samtools coverage $trimbam > ${trimbam.simpleName}.coverage.txt
	"""

}

process calculateRxy {

	// Calculate Rxy statistics for each individual
	
	publishDir "$params.outdir/06_Rxy", mode: 'copy'
	
	input:
	path(trimbam)
	path(rx_script)
	
	output:
	path("${trimbam.simpleName}.mapQ30.bam")
	path("${trimbam.simpleName}.mapQ30.coverage.txt")
	path("${trimbam.simpleName}.Rxy.txt")
	
	script:
	samtools_extra_threads = task.cpus - 1
	"""
	samtools view -q 30 -@ ${samtools_extra_threads} -o ${trimbam.simpleName}.mapQ30.bam $trimbam
	samtools coverage ${trimbam.simpleName}.mapQ30.bam > ${trimbam.simpleName}.mapQ30.coverage.txt
	Rscript $rx_script ${trimbam.simpleName}.mapQ30.coverage.txt > ${trimbam.simpleName}.Rxy.txt
	"""
	
}

process kmerSex {

	// Estimate sex assignment using kmers
	publishDir "$params.outdir/07_KmerSex", mode: 'copy'
	
	input:
	tuple val(sample), val(library), path(reads), val(rg)
	path(kmers)
	path(refseq)
	path "*"
	val(sry)
	
	output:
	path("${sample}_kmer.fq.gz")
	path("${sample}_kmer.bam")
	path("${sample}_kmer.cov_out")
	path("${sample}_kmer.sdry.cov")
	
	"""
	if [ -f *truncated.gz ]; then for i in *truncated.gz; do ln -s \${i} \${i%.gz}.fastq.gz; done; fi
	cat *fastq.gz > ${sample}.fq.gz
	bbduk.sh in=${sample}.fq.gz ref=${kmers} outm=${sample}_kmer.fq k=21
	bwa mem -M ${refseq} ${sample}_kmer.fq | samtools sort - | samtools markdup - ${sample}_kmer.bam
	samtools index ${sample}_kmer.bam
	samtools coverage -o ${sample}_kmer.cov_out -m -w 100 ${sample}_kmer.bam
	samtools coverage -r ${sry} -o ${sample}_kmer.sdry.cov ${sample}_kmer.bam 
	gzip ${sample}_kmer.fq
	"""
}

process blastUnalignedReads {

	// Extract unaligned reads in FASTA format for BLAST analysis
	// Blast unaligned reads against NT

	publishDir "$params.outdir/08_BLASTMetagenome", mode: 'copy'

	input:
	path finalbam
	path blastdb

	output:
	path "${finalbam.simpleName}.fa.gz"
	path "${finalbam.simpleName}.xml.gz"
	path "${finalbam.simpleName}.rma6"

	script:
	samtools_extra_threads = task.cpus - 1
	"""
	samtools fasta -@ ${samtools_extra_threads} -f 4 ${finalbam} > ${finalbam.simpleName}.fa
	blastn -query ${finalbam.simpleName}.fa -num_threads ${task.cpus} -db $blastdb -outfmt 5 -out >(gzip > ${finalbam.simpleName}.xml.gz)
	gzip ${finalbam.simpleName}.fa
	blast2rma -i ${finalbam.simpleName}.xml.gz -f BlastXML -bm BlastN -r ${finalbam.simpleName}.fa.gz -o ${finalbam.simpleName}.rma6
	"""

}

workflow {
	main:
		prepareRef(params.refseq)
		genMapIndex(params.refseq, params.gm_tmpdir) | genMapMap
		pe_read_data = Channel.fromPath(params.pelibraries).splitCsv(header:true).map { row -> tuple(row.Sample, row.Library, file(params.readDir + row.Read1), file(params.readDir + row.Read2), row.Adapter1, row.Adapter2, '@RG\\tID:' + row.Library + '\\tSM:' + row.Sample + '\\tLB:ILLUMINA\\tPL:ILLUMINA') }
		se_read_data = Channel.fromPath(params.selibraries).splitCsv(header:true).map { row -> tuple(row.Sample, row.Library, file(params.readDir + row.Read1), row.Adapter1, row.Adapter2, '@RG\\tID:' + row.Library + '\\tSM:' + row.Sample + '\\tLB:ILLUMINA\\tPL:ILLUMINA') }
		trimPEAdapters(pe_read_data)
		trimSEAdapters(se_read_data)
		all_reads = trimPEAdapters.out.mix(trimSEAdapters.out)
		alignSeqs(all_reads, params.refseq, prepareRef.out)
		realignIndels(alignSeqs.out.bam, alignSeqs.out.sample, params.refseq, prepareRef.out) | markDuplicates | profileDamage
		mergeLibraries(markDuplicates.out.groupTuple(by: 1)) // Need unique samples matched with their file paths
		reRealignIndels(mergeLibraries.out, params.refseq, prepareRef.out) | reMarkDuplicates | trimAncientTermini | calculateStatistics
		calculateRxy(trimAncientTermini.out, params.rx_script)
		kmerSex(all_reads.groupTuple(by: 0), params.kmers, params.refseq, prepareRef.out, params.sry)
		blastUnalignedReads(reMarkDuplicates.out, params.blastdb)
}
