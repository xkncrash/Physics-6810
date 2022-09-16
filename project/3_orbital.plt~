#  file: derivative_test.plt 
#
#  Gnuplot plot file for derivative_test output
#  
#  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
#               Ralf Bundschuh  bundshcuh.2@osu.edu
# 		Kyle Neumann	neumann.110@osu.edu
# 
#  Revision history
#   2004-01-24  original version for 780.20 session 5
#   2004-01-16  added postscript enhanced and comments for session 4
#   2021-12-29  made terminal type agnostic and restored output at the end
#   2022-01-26  modified version to give 3 fit functions

# record the time and date the graph was generated
set timestamp

# titles and labels
set title 'Test of Numerical Derivatives using alpha*x^beta'
set xlabel 'x'
set ylabel 'y'

# move the legend to a free space
set key left
set size ratio 1

# set the x and y axis scales (already logs)
set xrange [-2:2]
set yrange [-2:2]


# plot the data as well as the fit, with appropriate titles 
plot "2_orbital.dat" using 2:3 title "m1" with line, \
"2_orbital.dat" using 4:5 title "m2" with line

# output the plot to the file derivative_test.ps
# remember terminal type
set term push
set term postscript enhanced color
set out "2_orbital.ps"
replot
set out
# restore terminal type
set term pop
