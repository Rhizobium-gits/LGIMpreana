set terminal svg size 1600,900 enhanced font 'Arial,12' background '#FFFFFF'
set output '/tmp/microbiome_results/Figure11_DominantTaxa.svg'

set multiplot layout 3,3 title 'Dominant Taxa Dynamics Over Time' font 'Arial Bold,16'

set title 'Donor 1 / 0g' font 'Arial Bold,11'
set xlabel 'Time (hours)' font 'Arial,9'
set ylabel 'Abundance' font 'Arial,9'
set key top right font 'Arial,8'
set grid
set xrange [0:30]
set yrange [0:*]

$obs_D1_0g_0 << EOD
8.0 0.492635
16.0 0.463375
24.0 0.329807
EOD

$obs_D1_0g_1 << EOD
8.0 0.093993
16.0 0.126300
24.0 0.251178
EOD

$obs_D1_0g_2 << EOD
8.0 0.073844
16.0 0.063280
24.0 0.035912
EOD

plot $obs_D1_0g_0 using 1:2 title 'Bacteroides' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#E64B35', \
     $obs_D1_0g_1 using 1:2 title 'Pseudomonas' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#4DBBD5', \
     $obs_D1_0g_2 using 1:2 title 'Blautia' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#00A087'

set title 'Donor 1 / 1g' font 'Arial Bold,11'
set xlabel 'Time (hours)' font 'Arial,9'
set ylabel 'Abundance' font 'Arial,9'
set key top right font 'Arial,8'
set grid
set xrange [0:30]
set yrange [0:*]

$obs_D1_1g_0 << EOD
8.0 0.488325
16.0 0.504312
24.0 0.390492
EOD

$obs_D1_1g_1 << EOD
8.0 0.051920
16.0 0.115050
24.0 0.209161
EOD

$obs_D1_1g_2 << EOD
8.0 0.064509
16.0 0.048398
24.0 0.042393
EOD

plot $obs_D1_1g_0 using 1:2 title 'Bacteroides' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#E64B35', \
     $obs_D1_1g_1 using 1:2 title 'Pseudomonas' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#4DBBD5', \
     $obs_D1_1g_2 using 1:2 title 'Blautia' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#00A087'

set title 'Donor 1 / 5g' font 'Arial Bold,11'
set xlabel 'Time (hours)' font 'Arial,9'
set ylabel 'Abundance' font 'Arial,9'
set key top right font 'Arial,8'
set grid
set xrange [0:30]
set yrange [0:*]

$obs_D1_5g_0 << EOD
8.0 0.551234
16.0 0.571866
24.0 0.508788
EOD

$obs_D1_5g_1 << EOD
8.0 0.039920
16.0 0.071947
24.0 0.162339
EOD

$obs_D1_5g_2 << EOD
8.0 0.069826
16.0 0.053276
24.0 0.034336
EOD

plot $obs_D1_5g_0 using 1:2 title 'Bacteroides' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#E64B35', \
     $obs_D1_5g_1 using 1:2 title 'Pseudomonas' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#4DBBD5', \
     $obs_D1_5g_2 using 1:2 title 'Blautia' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#00A087'

set title 'Donor 2 / 0g' font 'Arial Bold,11'
set xlabel 'Time (hours)' font 'Arial,9'
set ylabel 'Abundance' font 'Arial,9'
set key top right font 'Arial,8'
set grid
set xrange [0:30]
set yrange [0:*]

$obs_D2_0g_0 << EOD
8.0 0.115494
16.0 0.108714
24.0 0.077666
EOD

$obs_D2_0g_1 << EOD
8.0 0.001034
16.0 0.000000
24.0 0.000050
EOD

$obs_D2_0g_2 << EOD
8.0 0.053308
16.0 0.044772
24.0 0.028552
EOD

plot $obs_D2_0g_0 using 1:2 title 'Bacteroides' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#E64B35', \
     $obs_D2_0g_1 using 1:2 title 'Pseudomonas' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#4DBBD5', \
     $obs_D2_0g_2 using 1:2 title 'Blautia' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#00A087'

set title 'Donor 2 / 1g' font 'Arial Bold,11'
set xlabel 'Time (hours)' font 'Arial,9'
set ylabel 'Abundance' font 'Arial,9'
set key top right font 'Arial,8'
set grid
set xrange [0:30]
set yrange [0:*]

$obs_D2_1g_0 << EOD
8.0 0.117384
16.0 0.135423
24.0 0.103078
EOD

$obs_D2_1g_1 << EOD
8.0 0.000000
16.0 0.000000
24.0 0.000068
EOD

$obs_D2_1g_2 << EOD
8.0 0.060462
16.0 0.052338
24.0 0.035742
EOD

plot $obs_D2_1g_0 using 1:2 title 'Bacteroides' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#E64B35', \
     $obs_D2_1g_1 using 1:2 title 'Pseudomonas' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#4DBBD5', \
     $obs_D2_1g_2 using 1:2 title 'Blautia' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#00A087'

set title 'Donor 2 / 5g' font 'Arial Bold,11'
set xlabel 'Time (hours)' font 'Arial,9'
set ylabel 'Abundance' font 'Arial,9'
set key top right font 'Arial,8'
set grid
set xrange [0:30]
set yrange [0:*]

$obs_D2_5g_0 << EOD
8.0 0.149599
16.0 0.128171
24.0 0.116809
EOD

$obs_D2_5g_1 << EOD
8.0 0.000029
16.0 0.000000
24.0 0.000000
EOD

$obs_D2_5g_2 << EOD
8.0 0.056597
16.0 0.050675
24.0 0.039647
EOD

plot $obs_D2_5g_0 using 1:2 title 'Bacteroides' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#E64B35', \
     $obs_D2_5g_1 using 1:2 title 'Pseudomonas' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#4DBBD5', \
     $obs_D2_5g_2 using 1:2 title 'Blautia' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#00A087'

set title 'Donor 3 / 0g' font 'Arial Bold,11'
set xlabel 'Time (hours)' font 'Arial,9'
set ylabel 'Abundance' font 'Arial,9'
set key top right font 'Arial,8'
set grid
set xrange [0:30]
set yrange [0:*]

$obs_D3_0g_0 << EOD
8.0 0.172831
16.0 0.181452
24.0 0.153220
EOD

$obs_D3_0g_1 << EOD
8.0 0.009628
16.0 0.014352
24.0 0.026075
EOD

$obs_D3_0g_2 << EOD
8.0 0.098767
16.0 0.096112
24.0 0.062731
EOD

plot $obs_D3_0g_0 using 1:2 title 'Bacteroides' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#E64B35', \
     $obs_D3_0g_1 using 1:2 title 'Pseudomonas' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#4DBBD5', \
     $obs_D3_0g_2 using 1:2 title 'Blautia' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#00A087'

set title 'Donor 3 / 1g' font 'Arial Bold,11'
set xlabel 'Time (hours)' font 'Arial,9'
set ylabel 'Abundance' font 'Arial,9'
set key top right font 'Arial,8'
set grid
set xrange [0:30]
set yrange [0:*]

$obs_D3_1g_0 << EOD
8.0 0.204062
16.0 0.241347
24.0 0.207330
EOD

$obs_D3_1g_1 << EOD
8.0 0.008795
16.0 0.014402
24.0 0.016600
EOD

$obs_D3_1g_2 << EOD
8.0 0.110155
16.0 0.077076
24.0 0.063140
EOD

plot $obs_D3_1g_0 using 1:2 title 'Bacteroides' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#E64B35', \
     $obs_D3_1g_1 using 1:2 title 'Pseudomonas' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#4DBBD5', \
     $obs_D3_1g_2 using 1:2 title 'Blautia' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#00A087'

set title 'Donor 3 / 5g' font 'Arial Bold,11'
set xlabel 'Time (hours)' font 'Arial,9'
set ylabel 'Abundance' font 'Arial,9'
set key top right font 'Arial,8'
set grid
set xrange [0:30]
set yrange [0:*]

$obs_D3_5g_0 << EOD
8.0 0.207659
16.0 0.247407
24.0 0.218245
EOD

$obs_D3_5g_1 << EOD
8.0 0.004378
16.0 0.006019
24.0 0.014884
EOD

$obs_D3_5g_2 << EOD
8.0 0.087338
16.0 0.100761
24.0 0.082786
EOD

plot $obs_D3_5g_0 using 1:2 title 'Bacteroides' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#E64B35', \
     $obs_D3_5g_1 using 1:2 title 'Pseudomonas' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#4DBBD5', \
     $obs_D3_5g_2 using 1:2 title 'Blautia' with linespoints pt 7 ps 1.2 lw 2 lc rgb '#00A087'

unset multiplot
