#! /bin/bash

echo "### 1D Shockwave Test ###"

# 1D Shockwave files
result=0
filenotfound=0
temp=0
fileamnt=$(find *1d*till* 2>/dev/null | wc -l )

if [ $fileamnt -ne 0 ]
then
	for file in *1d*till5*
	do
		diff -q $file ./DiffableResults/1DShock/$file &>/dev/null
		temp=$?
		if [ $temp -eq 2 ]
		then
			filenotfound=$(($filenotfound + 1))
		fi
		result=$(($result + $temp))
	done
fi

if [ $filenotfound -ne 0 ]
then
	if [ $fileamnt -ne $filenotfound ]
	then
		foundfiles=$(($fileamnt - $filenotfound))
		samefiles=$(($foundfiles - $(($result - ($filenotfound * 2)))))
		difffiles=$(($foundfiles - $samefiles))
		echo "$filenotfound/$fileamnt uncomparable 1DShock files"
		echo "$samefiles/$fileamnt 1DShock files are the same"
		echo "$difffiles/$fileamnt 1DShock files are different"
	else
		echo "No comparable 1DShock files not found"
	fi
elif [ $fileamnt -eq 0 ]
then
	echo "No 1DShock files found"
else
	if [ $result -eq 0 ]
	then
		echo -e "\xE2\x9C\x93 1DShock files are the same"
	else
		echo -e "\xE2\x9C\x97 1DShock files are different"
	fi
fi

echo ""
