set terminal svg size 1200,400 enhanced font 'Arial,12' background '#FFFFFF'
set output '/tmp/microbiome_results/Figure7_LinearPrediction.svg'

set multiplot layout 1,3 title 'Predicted 48h Composition by Gravity' font 'Arial Bold,14'

set title 'Microgravity (0g)' font 'Arial Bold,12'
set style fill solid 0.7
set boxwidth 0.7
set ylabel 'Predicted Abundance'
set xtics rotate by -45 right font 'Arial,9'
set grid y
set bars 2

$pred_0g << EOD
0 0.3383 0.1500 "Staphylo.."
1 0.1793 0.2193 "Pseudomo.."
2 0.1477 0.0954 "Escheric.."
3 0.0945 0.0873 "Acinetob.."
4 0.0784 0.0724 "Bacteroi.."
EOD
plot $pred_0g using 1:2:3:xtic(4) with boxerrorbars lc rgb '#E64B35' notitle

set title 'Earth (1g)' font 'Arial Bold,12'
set style fill solid 0.7
set boxwidth 0.7
set ylabel 'Predicted Abundance'
set xtics rotate by -45 right font 'Arial,9'
set grid y
set bars 2

$pred_1g << EOD
0 0.3176 0.1166 "Staphylo.."
1 0.1793 0.0984 "Bacteroi.."
2 0.1578 0.2052 "Pseudomo.."
3 0.0964 0.0707 "Escheric.."
4 0.0916 0.0987 "Acinetob.."
EOD
plot $pred_1g using 1:2:3:xtic(4) with boxerrorbars lc rgb '#00A087' notitle

set title 'Hypergravity (5g)' font 'Arial Bold,12'
set style fill solid 0.7
set boxwidth 0.7
set ylabel 'Predicted Abundance'
set xtics rotate by -45 right font 'Arial,9'
set grid y
set bars 2

$pred_5g << EOD
0 0.2844 0.1216 "Staphylo.."
1 0.2493 0.1641 "Bacteroi.."
2 0.1255 0.1584 "Pseudomo.."
3 0.0775 0.0635 "Escheric.."
4 0.0654 0.0643 "Acinetob.."
EOD
plot $pred_5g using 1:2:3:xtic(4) with boxerrorbars lc rgb '#F39B7F' notitle

unset multiplot
