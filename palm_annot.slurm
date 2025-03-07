#!/bin/bash -l
#SBATCH -o palm_annot.log
#SBATCH --job-name="palm_annot"
#SBATCH --nodes="1"
#SBATCH --ntasks-per-node="72"
#SBATCH --ntasks="72"
#SBATCH --time="12:00:00"
#SBATCH --account="lp_jm_virome_group"
#SBATCH -M wice

#You can download palm_annot.py/palmscan by cloning the git repo from https://github.com/rcedgar/palm_annot
#Add the bin and py directories to your PATH

cd ~/staging/Ludovica/phylogenetics/ictv_phylogenies

conda activate palm_annot

for i in *; do
	if [[ -d "$i" && "$i" != "orthomyxoviridae" ]]; then

		current_dir=$(pwd)

		cd $i || exit

		#Run palm_annot.py to identify viral RdRps (minscore of 90 == very high confidence)
		#I had to install a conda environment for 'hmmer' to have the right glibc version
		palm_annot.py --input $i.faa --seqtype aa --minscore 90 --rdrp $i.rdrp.faa --threads 72

		palmscan2 -search_pssms $i.rdrp.faa -core $i.core.faa -minflanklen 100 -threads 72

		cd "$current_dir" || exit
	fi
done
