#! /bin/bash

echo "### 2D Test ###"

# 2D files
result=0
filenotfound=0
temp=0
fileamnt=$(find *2d*till* 2>/dev/null | wc -l )

if [ $fileamnt -ne 0 ]
then
	for file in *2d*till5*
	do
		diff -q $file ./DiffableResults/$file 
		temp=$?
		if [ $temp -eq 2 ]
		then
			filenotfound=$(($filenotfound + 1))
		fi
		result=$(($result + $temp))
	done
fi

echo "$fileamnt 2D files compared"

if [ $filenotfound -ne 0 ]
then
	if [ $fileamnt -ne $filenotfound ]
	then
		foundfiles=$(($fileamnt - $filenotfound))
		samefiles=$(($foundfiles - $(($result - ($filenotfound * 2)))))
		difffiles=$(($foundfiles - $samefiles))
		echo "$filenotfound/$fileamnt uncomparable 2D files"
		echo "$samefiles/$fileamnt 2D files are the same"
		echo "$difffiles/$fileamnt 2D files are different"
	else
		echo "No comparable 2D files not found"
	fi
elif [ $fileamnt -eq 0 ]
then
	echo "No 2D Files found"
else
	if [ $result -eq 0 ]
	then
		echo -e "\xE2\x9C\x93 2D Files are the same"
	else
		echo -e "\xE2\x9C\x97 2D Files are different"
	fi
fi

echo ""
