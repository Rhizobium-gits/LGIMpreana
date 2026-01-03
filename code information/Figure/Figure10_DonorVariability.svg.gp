set terminal svg size 800,500 enhanced font 'Arial,12' background '#FFFFFF'
set output '/tmp/microbiome_results/Figure10_DonorVariability.svg'

set title 'Donor Variability Index Across Gravity Conditions' font 'Arial Bold,14'
set palette defined (0 '#F7FCF5', 0.2 '#C7E9C0', 0.4 '#74C476', 0.6 '#31A354', 0.8 '#006D2C', 1 '#00441B')
set cbrange [0:0.5]
set cblabel 'Variability Index' font 'Arial,10'
set xrange [-0.5:4.5]
set yrange [-0.5:2.5]

set xtics ('0g' 0, '1/6g' 1, '1g' 2, '1g-S' 3, '5g' 4) font 'Arial,10'
set ytics ('Donor 1' 0, 'Donor 2' 1, 'Donor 3' 2) font 'Arial,10'
set xlabel 'Gravity Condition' font 'Arial Bold,11'
set ylabel 'Donor' font 'Arial Bold,11'

$heatdata << EOD
0 0 -0.5 0.5 -0.5 0.5 0.1883
1 0 0.5 1.5 -0.5 0.5 0.1696
2 0 1.5 2.5 -0.5 0.5 0.1770
3 0 2.5 3.5 -0.5 0.5 0.2172
4 0 3.5 4.5 -0.5 0.5 0.1443
0 1 -0.5 0.5 0.5 1.5 0.1878
1 1 0.5 1.5 0.5 1.5 0.1533
2 1 1.5 2.5 0.5 1.5 0.1705
3 1 2.5 3.5 0.5 1.5 0.1971
4 1 3.5 4.5 0.5 1.5 0.1622
0 2 -0.5 0.5 1.5 2.5 0.2129
1 2 0.5 1.5 1.5 2.5 0.2257
2 2 1.5 2.5 1.5 2.5 0.2346
3 2 2.5 3.5 1.5 2.5 0.2229
4 2 3.5 4.5 1.5 2.5 0.2287
EOD

set style fill solid 1.0 noborder
plot $heatdata using 1:2:3:4:5:6:7 with boxxyerror palette notitle
