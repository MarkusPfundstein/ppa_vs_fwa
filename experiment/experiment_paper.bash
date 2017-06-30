# PURPOSE
# Compare ppalevy , ppa & fwa on standard functions with a lot of dimensions . 


NAME="experiment_paper"
ALGS="ppalevy ppa fwa"
FUNS="sphere rastrigrin ellipse cigar tablet schwefelFWA michalewicz12 ackley schwefel rosenbrock schwefel7 griewank"
DIMS="2 4 8 16 32 64 128"
SHIFTS="0 12 25"

LB=-100
UP=100
MAXG=8000
MAXFEVALS=600000
INITIALSIZE=10
NRUNS=50

echo "NRUNS: "$NRUNS
echo "FUNS: "$FUNS
echo "ALGS: "$ALGS
echo "DIMS: "$DIMS
echo "MAXG: "$MAXG
echo "MAXFEVLS: "$MAXFEVALS
echo "SHIFTS: "$SHIFTS

if [ ! -d $NAME ]; then 
	mkdir $NAME
fi

for d in $DIMS; do
	for a in $ALGS; do
		for f in $FUNS; do 
		    for s in $SHIFTS; do
				FEVALS=$MAXFEVALS
				if [ "$a" == "ppa" ] || [ "$a" == "ppalevy" ]; then
					echo "fix evals"
					FEVALS=$(($FEVALS * 3))
				fi
				echo "$a - $f - $d - $FEVALS"
				./solver-sop.exe \
					--origin-shift=$s \
					--algorithm=$a \
					--function=$f \
					--dimensions=$d \
					--ppa-nmax=50 \
					--max-generations=$MAXG \
					--max-fevals=$FEVALS \
					--min-bound=-100 \
					--max-bound=100 \
					--initial-size=$INITIALSIZE \
					--number-runs=$NRUNS \
					--write-values=$NAME"/"$a"_"$f"_"$d"_"$s"_"$MAXFEVALS".json"
			done
		done
	done
done
