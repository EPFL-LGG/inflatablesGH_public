for threshold in 0.800 0.775 0.750 0.725 0.700 0.675 0.650 0.625 0.600 0.575 0.550 0.525 0.500 0.475 0.450; do
    echo $threshold
    python sphere_cap.py $threshold > stdout_$threshold.txt 2> stderr_$threshold.txt
done
