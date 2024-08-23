#!/bin/bash -l
#SBATCH -A lp_jm_virome_group
#SBATCH -M wice
#SBATCH -N 1
#SBATCH -c 96
#SBATCH -t 3-00:00:00
#SBATCH -o blackflies_trees.log
#SBATCH --partition=batch_sapphirerapids

cd ~/staging/Ludovica/phylogenetics/new_phylogenies
conda activate phylogenetics

for i in *.dedupe.fasta; do
	sed -i 's/\*//g' $i

	mkdir -p ${i%.dedupe.fasta}
	cd ${i%.dedupe.fasta}

	if [[ $i != "genomoviridae_complete.dedupe.fasta" ]]; then
		palmscan2 -search_pssms ../$i -core ${i%.dedupe.fasta}.core.fasta
	else
		cp ../$i ${i%.dedupe.fasta}.core.fasta
	fi

	# einsi --thread 5 ${i%.dedupe.fasta}.core.fasta > ${i%.dedupe.fasta}.aligned.fasta
	# trimal -in ${i%.dedupe.fasta}.aligned.fasta -out ${i%.dedupe.fasta}.trimmed.fasta -gappyout -keepheader
	muscle -align ${i%.dedupe.fasta}.core.fasta -diversified -threads 96 -output ${i%.dedupe.fasta}.efa

	muscle -resample ${i%.dedupe.fasta}.efa -minconf 0 -output ${i%.dedupe.fasta}.resampled@

	#muscle -disperse ${i%.dedupe.fasta}.efa
	muscle -maxcc ${i%.dedupe.fasta}.efa -output ${i%.dedupe.fasta}.maxcc.afa

	muscle -strip_gappy_cols ${i%.dedupe.fasta}.maxcc.afa -output ${i%.dedupe.fasta}.maxcc.stripped.afa

	#iqtree -T AUTO -s ${i%.dedupe.fasta}.maxcc.stripped.afa --prefix ${i%.dedupe.fasta}.maxcc
	fasttree -lg -gamma -seed 1 -nosupport ${i%.dedupe.fasta}.maxcc.stripped.afa >  ${i%.dedupe.fasta}.maxcc.newick

	for iteration in {1..100}; do
		#iqtree -T AUTO -s ${i%.dedupe.fasta}.resampled$iteration --prefix ${i%.dedupe.fasta}.resampled$iteration
		fasttree -lg -gamma -seed 1 -nosupport ${i%.dedupe.fasta}.resampled$iteration > ${i%.dedupe.fasta}.resampled$iteration.newick
		cat ${i%.dedupe.fasta}.resampled$iteration.newick >> ${i%.dedupe.fasta}_resampled_trees.newick
	done

	newick -conf ${i%.dedupe.fasta}.maxcc.newick -trees ${i%.dedupe.fasta}_resampled_trees.newick -output ${i%.dedupe.fasta}_conftree.newick
	cd ~/staging/Ludovica/phylogenetics/new_phylogenies
done
