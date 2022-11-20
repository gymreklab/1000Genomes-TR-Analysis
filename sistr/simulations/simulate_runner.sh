for i in {11..20}
do
	qsub simulate_AFR.pbs -v a=2,b=$i
done


for i in {5..13}
do
	qsub simulate_AFR.pbs -v a=3,b=$i
done


for i in {7..10}
do
	qsub simulate_AFR.pbs -v a=4,b=$i
done
