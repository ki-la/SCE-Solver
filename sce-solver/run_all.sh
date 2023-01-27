#!/bin/sh
shopt -s nullglob

cargo build --release

configs=( 1 2 3 4 )
for CON in ${configs[@]}
do
    echo $CON
    for DIR in ../test_instances/final/*/;  # loop over all directories
    do
        echo $DIR
        for FILE in $DIR*.dimacs # loop over all .dimacs files in DIR
        do
            echo $FILE
            if [[ $CON == 1 ]]
            then
		        ./target/release/sce-solver $FILE 3600 0 2 0 1 1 1 -o final/final.csv  # execute BB solver w/ best config
            elif [[ $CON == 2 ]]
            then
                ./target/release/sce-solver $FILE 3600 0 2 0 1 1 -1 -o final/final_no_dr.csv  # execute BB solver w/o DR
		    elif [[ $CON == 3 ]]
		    then
			    timeout 3600 ./target/release/sce-solver $FILE 3600 1 2 0 0 0 1 -o final/final_ilp.csv  # execute ILP solver w/ DR
		    else
			    timeout 3600 ./target/release/sce-solver $FILE 3600 1 2 0 0 0 0 -o final/final_ilp_no_dr.csv  # execute ILP solver w/o DR
		    fi
        done        
    done
done

