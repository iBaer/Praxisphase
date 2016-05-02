for file in d_*till5*Steps
do
	gnuplot <<- EndExec
	set term png
	set zrange [700:820]
	set output "$file.png"
	splot "$file"
	exit
	EndExec
done

for file in p_*till5*Steps
do
	gnuplot <<- EndExec
	set term png
	set zrange [4.2e+08:5.2e+08]
	set output "$file.png"
	splot "$file"
	exit
	EndExec
done

for file in ux_*till5*Steps
do
	gnuplot <<- EndExec
	set term png
	set zrange [-110:110]
	set output "$file.png"
	splot "$file"
	exit
	EndExec
done

for file in uxr_*till5*Steps
do
	gnuplot <<- EndExec
	set term png
	set zrange [-100:100]
	set output "$file.png"
	splot "$file"
	exit
	EndExec
done
