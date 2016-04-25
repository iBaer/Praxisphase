#! /bin/bash

clear
echo "##### MULTIFLOW DIFF TEST #####"

rm -f *till5*

loop=0
for(( fluxmethod=1 ; fluxmethod<=2 ; fluxmethod++ ))
do
	for(( initmethod=0 ; initmethod<=3 ; initmethod++ ))
	do
		for(( splitmethod=1 ; splitmethod<=2 ; splitmethod++ ))
		do
			./main > /dev/null <<- EndExec
			$fluxmethod
			$initmethod
			$splitmethod
			1
			EndExec
			loop=$(($loop + 1))
		done
	done
done

./diffable2D.sh

rm *2d*till5*

mv ./gas-liquid.in ./gas-liquid-2d.in
mv ./gas-liquid-1d.in ./gas-liquid.in

for(( initmethod=0 ; initmethod<=1 ; initmethod++ ))
do
	for(( fluxmethod=1 ; fluxmethod<=2 ; fluxmethod++ ))
	do
		./main > /dev/null <<- EndExec
		$fluxmethod
		$initmethod
		1
		EndExec
		loop=$(($loop + 1))
	done
	if [ $initmethod -eq 0 ]
	then
		./diffable1Drare.sh
		rm *1d*till5*
	else
		./diffable1Dshock.sh
		rm *1d*till5*
	fi
done

mv ./gas-liquid.in ./gas-liquid-1d.in
mv ./gas-liquid-2d.in ./gas-liquid.in

echo "##### ALL TESTS COMPLETE #####"
