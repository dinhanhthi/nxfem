# load "ffpp-fitted-cr.gnu"
#set logscale x 10
#set logscale y 10
#plot "ffpp-data.dat" using 1:2 w lp t "L2", "ffpp-data.dat" using 1:3 w lp t "H1"
plot "ffpp-data.dat" using (log10($1)):(log10($2)) w lp t "L2", "ffpp-data.dat" using (log10($1)):(log10($3)) w lp t "H1"
f(x) = a*x+b
set fit quiet #don't print on screen, print to "fit.log"
fit f(x) "ffpp-data.dat" using (log10($1)):(log10($2)) via a,b
g(x) = c*x+d
set fit quiet #don't print on screen, print to "fit.log"
fit g(x) "ffpp-data.dat" using (log10($1)):(log10($3)) via c,d
