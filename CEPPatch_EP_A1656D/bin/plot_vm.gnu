set terminal pngcairo size 1920,1080 font ",30,bold" linewidth 2
set output "plot_result500.png"

set title "AP Graph A1656D for TN2004EPI Scale 500" font ",40"
set xlabel "Time (msec)"
set ylabel "AP (Vm)"


plot "vmcheck_wt.plt" using 1:2 with lines title "Wildtype", "vmcheck_mt.plt" using 1:2 with lines title "Mutation", "vmcheck_me.plt" using 1:2 with lines title "Mexilatine", "vmcheck_fl.plt" using 1:2 with lines title "Flecinadine", "vmcheck_ra.plt" using 1:2 with lines title "Ranolazine"
