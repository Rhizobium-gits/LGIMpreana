set terminal svg size 1200,800 enhanced font 'Arial,12' background '#FFFFFF'
set output '/tmp/microbiome_results/Figure12_NetworkStructure.svg'

set multiplot layout 2,2 title 'Network Structure Metrics by Gravity' font 'Arial Bold,16'
set style fill solid 0.8 border -1
set boxwidth 0.7
set grid y

$netdata << EOD
"0g" 0.3289 0.0000 0.0504 0.0000 2485.0 2211.0
"1/6g" 0.3289 0.0000 0.0504 0.0000 2485.0 2211.0
"1g" 0.3289 0.0000 0.0504 0.0000 2485.0 2211.0
"5g" 0.3289 0.0000 0.0504 0.0000 2485.0 2211.0
EOD

set title 'Network Connectance' font 'Arial Bold,12'
set ylabel 'Connectance (±SD)' font 'Arial,10'
set bars 2
plot $netdata using 0:2:3:xtic(1) with boxerrorbars lc rgb '#4DBBD5' notitle

set title 'Mean Interaction Strength' font 'Arial Bold,12'
set ylabel 'Strength (±SD)' font 'Arial,10'
plot $netdata using 0:4:5:xtic(1) with boxerrorbars lc rgb '#00A087' notitle

set title 'Positive Interactions' font 'Arial Bold,12'
set ylabel 'Count' font 'Arial,10'
unset bars
plot $netdata using 0:6:xtic(1) with boxes lc rgb '#3C5488' notitle

set title 'Negative Interactions' font 'Arial Bold,12'
set ylabel 'Count' font 'Arial,10'
plot $netdata using 0:7:xtic(1) with boxes lc rgb '#E64B35' notitle

unset multiplot
