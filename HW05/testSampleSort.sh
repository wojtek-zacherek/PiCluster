for n in {1..30};
do
	qsub -I -l nodes=1:ppn=$n
	
done
