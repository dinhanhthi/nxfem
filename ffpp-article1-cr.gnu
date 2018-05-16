# using by call: load "ffpp-article1-cr.gnu"
# set logscale x 10
# set logscale y 10
#plot "ffpp-article1-cr.dat" using 1:2 w lp t "L2"
set ylabel 'log(err)'
set xlabel 'log(h)'
plot "ffpp-article1-cr.dat" using (log10($1)):(log10($2)) w lp t "L2v", '' using (log10($1)):(log10($3)) w lp t "H1v"

f(x) = a*x+b
# set fit quiet #don't print on screen, print to "fit.log"
fit f(x) "ffpp-article1-cr.dat" using (log10($1)):(log10($5)) via a,b
