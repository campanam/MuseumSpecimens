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
	tuple val(library), val(sample), path("${library}.fastq.gz"), val(rg)
	
	"""
	AdapterRemoval --file1 $reads1 --file2 $reads2 --basename $library --adapter1 $adapter1 --adapter2 $adapter2 --gzip --collapse --minlength 30
	cat ${library}.collapsed.gz ${library}.collapsed.truncated.gz > ${library}.fastq.gz
	"""

}

process trimSEAdapters {

	// Trim adapters for paired-end libraries using AdapterRemoval 2.3.3

	input:
	tuple val(sample), val(library), path(reads1), val(adapter1), val(adapter2), val(rg)
	
	output:
	tuple val(library), val(sample), path("${library}.fastq.gz"), val(rg)
	
	"""
	AdapterRemoval --file1 $reads1 --basename $library --adapter1 $adapter1 --adapter2 $adapter2 --gzip --collapse --minlength 30
	cat ${library}.collapsed.gz ${library}.collapsed.truncated.gz > ${library}.fastq.gz
	"""

}


// START HERE with alignSeqs and markDuplicates


bwa samse -r RGtag refseq.fa <(bwa aln -t number seqs) seqs | something















workflow {
	main:
		prepareRef(params.refseq)
		genMapIndex(params.refseq, params.gm_tmpdir) | genMapMap
		pe_read_data = Channel.fromPath(params.pelibraries).splitCsv(header:true).map { row -> tuple(row.Sample, row.Library, file(params.readDir + row.Read1), file(params.readDir + row.Read2), row.Adapter1, row.Adapter2, '@RG\\tID:' + row.Library + '\\tSM:' + row.Sample + '\\tLB:ILLUMINA\\tPL:ILLUMINA') }
		se_read_data = Channel.fromPath(params.pelibraries).splitCsv(header:true).map { row -> tuple(row.Sample, row.Library, file(params.readDir + row.Read1), row.Adapter1, row.Adapter2, '@RG\\tID:' + row.Library + '\\tSM:' + row.Sample + '\\tLB:ILLUMINA\\tPL:ILLUMINA') }
		trimPEAdapters(pe_read_data)
		trimSEAdapters(se_read_data)
		all_reads = trimPEAdapters.out.mix(trimSEAdapters.out)
		alignSeqs(read_data, params.refseq, prepareRef.out) | markDuplicates

}