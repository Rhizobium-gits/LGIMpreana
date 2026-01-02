set terminal svg size 1000,600 enhanced font 'Arial,12' background '#FFFFFF'
set output '/tmp/microbiome_results/Figure10_DonorVariability.svg'

set title 'Donor Variability Index by Gravity' font 'Arial Bold,14'
set palette defined (0 '#FFFFCC', 0.25 '#FEB24C', 0.5 '#FC4E2A', 1 '#BD0026')
set cbrange [0:0.5]
set cblabel 'Variability Index' font 'Arial,11'

$matrix << EOD
0.1170 0.1461 0.1287 0.1560 0.1093 
0.1452 0.1621 0.1371 0.1450 0.1380 
0.1840 0.1845 0.1727 0.1803 0.1879 
EOD

set xtics ('Microgravity' 0, 'Lunar (1/6g)' 1, 'Earth (1g)' 2, 'Static' 3, 'Hypergravity' 4) font 'Arial,10'
set ytics ('D1' 0, 'D2' 1, 'D3' 2) font 'Arial,10'
set xlabel 'Gravity Condition' font 'Arial Bold,12'
set ylabel 'Donor' font 'Arial Bold,12'

plot $matrix matrix with image notitle
