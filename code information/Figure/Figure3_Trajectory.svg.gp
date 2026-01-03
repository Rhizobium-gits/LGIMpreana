set terminal svg size 1600,500 enhanced font 'Arial,12' background '#FFFFFF'
set output '/tmp/microbiome_results/Figure3_Trajectory.svg'

set multiplot layout 1,3 title 'Temporal Trajectory in PCoA Space (8h→16h→24h)' font 'Arial Bold,16'

set title 'Donor 1' font 'Arial Bold,14'
set xlabel 'PC1' font 'Arial,11'
set ylabel 'PC2' font 'Arial,11'
set key outside right top font 'Arial,9' box
set grid

$traj_D1_0g << EOD
-0.297412 0.000000
-0.290323 0.000000
-0.244102 0.000000
EOD

$traj_D1_1g << EOD
-0.300887 0.000000
-0.299317 0.000000
-0.265489 0.000000
EOD

$traj_D1_5g << EOD
-0.300615 0.000000
-0.300603 0.000000
-0.297232 0.000000
EOD

plot $traj_D1_0g using 1:2 title '0g' with linespoints pt 7 ps 2 lw 2.5 lc rgb '#E64B35', \
     $traj_D1_1g using 1:2 title '1g' with linespoints pt 7 ps 2 lw 2.5 lc rgb '#00A087', \
     $traj_D1_5g using 1:2 title '5g' with linespoints pt 7 ps 2 lw 2.5 lc rgb '#F39B7F'

set title 'Donor 2' font 'Arial Bold,14'
set xlabel 'PC1' font 'Arial,11'
set ylabel 'PC2' font 'Arial,11'
set key outside right top font 'Arial,9' box
set grid

$traj_D2_0g << EOD
0.360033 0.000000
0.361392 0.000000
0.421856 0.000000
EOD

$traj_D2_1g << EOD
0.355361 0.000000
0.353338 0.000000
0.399144 0.000000
EOD

$traj_D2_5g << EOD
0.345898 0.000000
0.355952 0.000000
0.373448 0.000000
EOD

plot $traj_D2_0g using 1:2 title '0g' with linespoints pt 7 ps 2 lw 2.5 lc rgb '#E64B35', \
     $traj_D2_1g using 1:2 title '1g' with linespoints pt 7 ps 2 lw 2.5 lc rgb '#00A087', \
     $traj_D2_5g using 1:2 title '5g' with linespoints pt 7 ps 2 lw 2.5 lc rgb '#F39B7F'

set title 'Donor 3' font 'Arial Bold,14'
set xlabel 'PC1' font 'Arial,11'
set ylabel 'PC2' font 'Arial,11'
set key outside right top font 'Arial,9' box
set grid

$traj_D3_0g << EOD
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
EOD

$traj_D3_1g << EOD
-0.079681 0.000000
-0.123200 0.000000
0.000000 0.000000
EOD

$traj_D3_5g << EOD
-0.084484 0.000000
-0.097833 0.000000
0.000000 0.000000
EOD

plot $traj_D3_0g using 1:2 title '0g' with linespoints pt 7 ps 2 lw 2.5 lc rgb '#E64B35', \
     $traj_D3_1g using 1:2 title '1g' with linespoints pt 7 ps 2 lw 2.5 lc rgb '#00A087', \
     $traj_D3_5g using 1:2 title '5g' with linespoints pt 7 ps 2 lw 2.5 lc rgb '#F39B7F'

unset multiplot
