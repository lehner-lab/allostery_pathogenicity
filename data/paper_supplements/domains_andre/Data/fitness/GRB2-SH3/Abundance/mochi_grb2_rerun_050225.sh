#BSUB -cwd /lustre/scratch126/gengen/projects/alpha-allostery-global/git-p1-allo-seq-evo/data/paper_supplements/domains_andre/Data/fitness/GRB2-SH3/Abundance
#BSUB -o ./mochiRun.%J.out
#BSUB -e ./mochiRun.%J.err
#BSUB -J grb2-rerun
#BSUB -G team354
#BSUB -n 10 	# number of CPUs. Default: 1 # it doesn't require a lot of resources
#BSUB -q "normal"
#BSUB -R "select[mem>200000] rusage[mem=200000] span[hosts=1]" 	# RAM memory part 1. Default: 100MB
#BSUB -M200000 # # RAM memory part 2. Default: 100MB
#BSUB -W720 # time for the job HH:MM. 500 min

module load HGI/softpack/groups/team354/mochi/0.9

run_mochi.py  --model_design=model_design.txt 
