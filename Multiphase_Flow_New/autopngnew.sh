for file in d_*till5*Steps
do
	gnuplot <<- EndExec
	set term png
	set zrange [0:100]
	set view 90,0
	set output "$file.png"
	splot "$file"
	exit
	EndExec
done
