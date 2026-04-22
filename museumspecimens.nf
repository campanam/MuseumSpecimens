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
	
	publishDir  "$params.outdir/00_TrimmedReads", mode: 'copy', enabled: "$params.keep_trimmed_reads"

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
	
	publishDir  "$params.outdir/00_TrimmedReads", mode: 'copy', enabled: "$params.keep_trimmed_reads"

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
	
	script:
	if ( params.csi )
		"""
		samtools index -c ${rg_bam}
		$gatk LeftAlignIndels -R ${refseq} -I $rg_bam -O ${rg_bam.simpleName}.realn.bam --create-output-bam-index false
		"""
	else
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

	// Profile deamination damage for ancient libraries and calculate individual library statistics
	
	publishDir "$params.outdir/02_LibraryStatistics_DamageProfiles", mode: 'copy', pattern: "*_damage/*.*"
	publishDir "$params.outdir/02_LibraryStatistics_DamageProfiles", mode: 'copy', pattern: "*.markdup.stats.txt"
	publishDir "$params.outdir/02_LibraryStatistics_DamageProfiles", mode: 'copy', pattern: "*.markdup.depth.txt.gz"
	publishDir "$params.outdir/02_LibraryStatistics_DamageProfiles", mode: 'copy', pattern: "*.markdup.coverage.txt"
	
	input:
	tuple path(mrkdupbam), val(sample)
	path refseq
	path "*"
	
	output:
	path "${mrkdupbam.simpleName}_damage/*pdf"
	path "${mrkdupbam.simpleName}_damage/*txt"
	path "${mrkdupbam.simpleName}_damage/*log"
	path "${trimbam.simpleName}.markdup.*.txt*"
	
	script:
	samtools_extra_threads = task.cpus - 1
	"""
	samtools flagstat -@ ${samtools_extra_threads} $trimbam > ${trimbam.simpleName}.markdup.stats.txt
	samtools depth -@ ${samtools_extra_threads} $trimbam | gzip > ${trimbam.simpleName}.markdup.depth.txt.gz
	samtools coverage $trimbam > ${trimbam.simpleName}.markdup.coverage.txt
	damageprofiler ${params.java11_options} -i ${mrkdupbam} -o ${mrkdupbam.simpleName}_damage -r ${refseq}
	"""

}

process mergeLibraries {

	// Merge libraries by their sample IDs using SAMtools merge
	
	input:
	tuple path(bam), val(sample)
	
	output:
	tuple path("${sample}_merged.bam"), val(bams)
	
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
	tuple path(rg_bam, bams)
	path refseq
	path "*"
	
	output:
	tuple path("${rg_bam.simpleName}.realn.bam"), val(bams)
	
	script:
	if (bams == 1) // Skip realignment
		"""
		ln -s $rg_bam ${rg_bam.simpleName}.realn.bam
		"""
	else if ( params.csi )
		"""
		samtools index -c ${rg_bam}
		$gatk LeftAlignIndels -R ${refseq} -I $rg_bam -O ${rg_bam.simpleName}.realn.bam --create-output-bam-index false
		"""
	else
		"""
		$picard BuildBamIndex I=${rg_bam}
		$gatk LeftAlignIndels -R ${refseq} -I $rg_bam -O ${rg_bam.simpleName}.realn.bam
		"""

}

process reMarkDuplicates {

	// Mark duplicates using sambamba, samtools or picard
	
	publishDir  "$params.outdir/03_MergedBAMs", mode: 'copy'
	
	input:
	tuple path(sorted_bam), val(bams)
	
	output:
	path("${sorted_bam.simpleName}.markdup.bam")
	
	script:
	samtools_extra_threads = task.cpus - 1
	if (bams == 1) // Skip realignment
		"""
		ln -s $sorted_bam ${sorted_bam.simpleName}.realn.bam
		"""
	else if ( params.markDuplicates == "sambamba" )
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
	
	publishDir "$params.outdir/05_MergedStats", mode: 'copy', pattern: "*.markdup.stats.txt"
	publishDir "$params.outdir/05_MergedStats", mode: 'copy', pattern: "*.markdup.depth.txt.gz"
	publishDir "$params.outdir/05_MergedStats", mode: 'copy', pattern: "*.markdup.coverage.txt"
	publishDir "$params.outdir/07_MapQStats", mode: 'copy', pattern: "*.mapq.stats.txt"
	publishDir "$params.outdir/07_MapQStats", mode: 'copy', pattern: "*.mapq.depth.txt.gz"
	publishDir "$params.outdir/07_MapQStats", mode: 'copy', pattern: "*.mapq.coverage.txt"
	
	input:
	path(trimbam)
	val extension
	
	output:
	path "${trimbam.simpleName}.${extension}.*.txt*"
	
	script:
	samtools_extra_threads = task.cpus - 1
	"""
	samtools flagstat -@ ${samtools_extra_threads} $trimbam > ${trimbam.simpleName}.${extension}.stats.txt
	samtools depth -@ ${samtools_extra_threads} $trimbam | gzip > ${trimbam.simpleName}.${extension}.depth.txt.gz
	samtools coverage $trimbam > ${trimbam.simpleName}.${extension}.coverage.txt
	"""

}

process filterMapQ {

	// Filter by minimum MapQ
	
	publishDir "$params.outdir/06_MapQBAMs", mode: 'copy'
	
	input:
	path mrkdupbam
	val mapq
	
	output:
	path "${mrkdupbam.simpleName}.mapq.bam"
	
	script:
	samtools_extra_threads = task.cpus - 1
	if (mapq < 1)
		"""
		ln -s ${mrkdupbam} ${mrkdupbam.simpleName}.mapq.bam
		"""
	else
		"""
		samtools view -q $mapq -b -o ${mrkdupbam.simpleName}.mapq.bam $mrkdupbam
		"""
	
}

process calculateRxy {

	// Calculate Rxy statistics for each individual
	
	publishDir "$params.outdir/08_Rxy", mode: 'copy'
	
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
	publishDir "$params.outdir/09_KmerSex", mode: 'copy'
	
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

process extractUnalignedReads {

	// Extract unaligned reads in FASTA format for BLAST analysis
	// Remove duplicates using CD-HIT-EST

	publishDir "$params.outdir/10_BLASTMetagenome", mode: 'copy'

	input:
	path finalbam

	output:
	path "${finalbam.simpleName}.fa.gz"

	script:
	samtools_extra_threads = task.cpus - 1
	"""
	samtools fasta -@ ${samtools_extra_threads} -f 4 ${finalbam} > ${finalbam.simpleName}.nonuniq.fa
	cd-hit-est -c 1 -M 0 -i ${finalbam.simpleName}.nonuniq.fa -T ${task.cpus} -o ${finalbam.simpleName}.fa
	rm ${finalbam.simpleName}.nonuniq.fa
	rm *clstr
	gzip ${finalbam.simpleName}.fa
	"""

}

process blastUnalignedReads {

	// Blast collapsed unaligned reads against NT

	publishDir "$params.outdir/10_BLASTMetagenome", mode: 'copy'

	input:
	path uniqreads
	val blastdb

	output:
	path "${uniqreads.simpleName}.xml.gz"
	path "${uniqreads.simpleName}.rma6"

	"""
	blastn -query <(gunzip -fc $uniqreads) -db $blastdb -outfmt 5 -num_threads ${task.cpus} -out >(gzip > ${uniqreads.simpleName}.xml.gz)
	blast2rma -i ${uniqreads.simpleName}.xml.gz -f BlastXML -bm BlastN -r ${uniqreads} -o ${uniqreads.simpleName}.rma6
	"""

}

workflow mapqStats {
	// Alias of calculateStatistics for MapQ-filtered data
	take:
		bam
	main:
		calculateStatistics(bam, "mapq")
	emit:
		calculateStatistics.out
}

workflow {
	main:
		prepareRef(params.refseq)
		if (params.genmap) { genMapIndex(params.refseq, params.gm_tmpdir) | genMapMap }
		if (params.pelibraries != "NULL") {
			pe_read_data = Channel.fromPath(params.pelibraries).splitCsv(header:true).map { row -> tuple(row.Sample, row.Library, file(params.readDir + row.Read1), file(params.readDir + row.Read2), row.Adapter1, row.Adapter2, '@RG\\tID:' + row.RG + '\\tSM:' + row.Sample + '\\tLB:' + row.Library + '\\tPL:' + row.PL) }
			trimPEAdapters(pe_read_data)
		}
		if (params.selibraries != "NULL") {
			se_read_data = Channel.fromPath(params.selibraries).splitCsv(header:true).map { row -> tuple(row.Sample, row.Library, file(params.readDir + row.Read1), row.Adapter1, row.Adapter2, '@RG\\tID:' + row.RG + '\\tSM:' + row.Sample + '\\tLB:' + row.Library + '\\tPL:' + row.PL') }
			trimSEAdapters(se_read_data)
		}
		if (params.pelibraries != "NULL" && params.selibraries != "NULL") {
			all_reads = trimPEAdapters.out.mix(trimSEAdapters.out)
		} else if (params.pelibraries == "NULL") {
			all_reads = trimSEAdapters.out
		} else {
			all_reads = trimPEAdapters.out
		}
		alignSeqs(all_reads, params.refseq, prepareRef.out)
		realignIndels(alignSeqs.out.bam, alignSeqs.out.sample, params.refseq, prepareRef.out) | markDuplicates 
		profileDamage(markDuplicates.out, params.refseq, prepareRef.out)
		mergeLibraries(markDuplicates.out.groupTuple(by: 1)) // Need unique samples matched with their file paths
		reRealignIndels(mergeLibraries.out, params.refseq, prepareRef.out) | reMarkDuplicates | trimAncientTermini 
		calculateStatistics(trimAncientTermini.out,"markdup")
		if (params.mapq > 0) {
			filterMapQ(trimAncientTermini.out, params.mapq)
			mapqStats(filterMapQ.out)
		}
		if (params.rx) { calculateRxy(filterMapQ.out, params.rx_script) }
		if (params.kmerSex) { kmerSex(all_reads.groupTuple(by: 0), params.kmers, params.refseq, prepareRef.out, params.sry) }
		if (params.blast) {
			extractUnalignedReads(reMarkDuplicates.out)
			blastUnalignedReads(extractUnalignedReads.out, params.blastdb)
		}
}
