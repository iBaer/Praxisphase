#! /bin/bash

clear
echo "##### MULTIFLOW DIFF TEST #####"

cd ./DiffableResults/ 
echo "Total comparable 2D files $(ls *till* | wc -l)"
cd ./1DRare/
echo "Total comparable 1DRare files $(ls *till* | wc -l)"
cd ..
cd ./1DShock/
echo "Total comparable 1DShock files $(ls *till* | wc -l)"
echo ""
cd ..
cd ..

rm -f *till5*

START=$(date +%s%N)
totalfiles=0

loop=0
for(( initmethod=0 ; initmethod<=4 ; initmethod++ ))
do
	for(( fluxmethod=0 ; fluxmethod<=1 ; fluxmethod++ ))
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
	for(( fluxmethod=0 ; fluxmethod<=1 ; fluxmethod++ ))
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

END=$(date +%s%N)
DIFF=$(( $END - $START ))
echo "Tests finished in $(( DIFF / 1000000000 ))s $(( DIFF / 1000000 ))ms"
mv ./gas-liquid.in ./gas-liquid-1d.in
mv ./gas-liquid-2d.in ./gas-liquid.in

echo "##### ALL TESTS COMPLETE #####"

./diffstart-analytisch.sh
