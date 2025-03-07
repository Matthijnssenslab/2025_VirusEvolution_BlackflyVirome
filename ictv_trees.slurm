#!/bin/bash -l
#SBATCH -A lp_jm_virome_group
#SBATCH -M wice
#SBATCH -N 1
#SBATCH -c 96
#SBATCH -t 3-00:00:00
#SBATCH -o ictv_trees_monodnaviria.log
#SBATCH --partition=batch_sapphirerapids

cd ~/staging/Ludovica/phylogenetics/ictv_phylogenies
conda activate phylogenetics

#for i in *; do
for i in parvoviridae genomoviridae; do
	if [ -d $i ]; then

	cd $i

	muscle -align ${i}.core.faa -diversified -threads 96 -output ${i}.efa

	muscle -resample ${i}.efa -minconf 0 -output ${i}.resampled@

	#muscle -disperse ${i}.efa
	muscle -maxcc ${i}.efa -output ${i}.maxcc.afa

	muscle -strip_gappy_cols ${i}.maxcc.afa -output ${i}.maxcc.stripped.afa

	#iqtree -T AUTO -s ${i}.maxcc.stripped.afa --prefix ${i}.maxcc
	fasttree -lg -gamma -seed 1 -nosupport ${i}.maxcc.stripped.afa >  ${i}.maxcc.newick

	for iteration in {1..100}; do
		#iqtree -T AUTO -s ${i}.resampled$iteration --prefix ${i}.resampled$iteration
		fasttree -lg -gamma -seed 1 -nosupport ${i}.resampled$iteration > ${i}.resampled$iteration.newick
		cat ${i}.resampled$iteration.newick >> ${i}_resampled_trees.newick
	done

	newick -conf ${i}.maxcc.newick -trees ${i}_resampled_trees.newick -output ${i}_conftree.newick

	cd ..

	fi
done
