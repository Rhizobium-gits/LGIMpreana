set terminal svg size 1000,700 enhanced font 'Arial,12' background '#FFFFFF'
set output '/tmp/microbiome_results/Figure12_NetworkStructure.svg'

set multiplot layout 2,2 title 'Network Structure Across Gravity' font 'Arial Bold,14'
set style fill solid 0.7
set boxwidth 0.8
set xtics rotate by -30 right font 'Arial,9'
set grid y

$data << EOD
"Microgravity" 0.0000 0.0000 0.0 0.0
"Lunar (1/6g)" 0.0000 0.0000 0.0 0.0
"Earth (1g)" 0.0000 0.0000 0.0 0.0
"Static" 0.0000 0.0000 0.0 0.0
"Hypergravity" 0.0000 0.0000 0.0 0.0
EOD

set title 'Connectance'
plot $data using 0:2:xtic(1) with boxes lc rgb '#4DBBD5' notitle

set title 'Mean Interaction Strength'
plot $data using 0:3:xtic(1) with boxes lc rgb '#00A087' notitle

set title 'Positive Interactions'
plot $data using 0:4:xtic(1) with boxes lc rgb '#3C5488' notitle

set title 'Negative Interactions'
plot $data using 0:5:xtic(1) with boxes lc rgb '#E64B35' notitle

unset multiplot
