set terminal svg size 900,700 enhanced font 'Arial,12' background '#FFFFFF'
set output '/tmp/microbiome_results/Figure3_Trajectory.svg'

set title 'Temporal Trajectory in PCoA Space' font 'Arial Bold,14'
set xlabel 'PC1' font 'Arial Bold,12'
set ylabel 'PC2' font 'Arial Bold,12'
set key outside right top box font 'Arial,10'
set grid

$traj_0g << EOD
0.020874 0.000000
0.023690 0.000000
0.059251 0.000000
EOD

$traj_1g << EOD
-0.008403 0.000000
-0.023059 0.000000
0.044552 0.000000
EOD

$traj_5g << EOD
-0.013067 0.000000
-0.014161 0.000000
0.025405 0.000000
EOD

plot $traj_0g using 1:2 title 'Microgravity' with linespoints pt 7 ps 1.5 lw 2 lc rgb '#E64B35', \
     $traj_1g using 1:2 title 'Earth (1g)' with linespoints pt 7 ps 1.5 lw 2 lc rgb '#00A087', \
     $traj_5g using 1:2 title 'Hypergravity' with linespoints pt 7 ps 1.5 lw 2 lc rgb '#F39B7F'
