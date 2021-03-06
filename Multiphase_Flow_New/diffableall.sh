#! /bin/bash

# 2D files
result=0
filenotfound=0
temp=0
for file in *2d*till5*
do
	diff -q $file ./DiffableResults/$file &>/dev/null
	temp=$?
	if [ $temp -eq 2 ]
	then
		filenotfound=$(($filenotfound + 1))
	fi
	result=$(($result + $temp))
done

if [ $filenotfound -ne 0 ]
then
	if [ $result -ne $(($filenotfound * 2)) ]
	then
		echo "$filenotfound 2D Files not found"
	else
		echo "2D Files not found"
	fi
else
	if [ $result -eq 0 ]
	then
		echo "2D Files are the same"
	else
		echo "2D Files are different"
	fi
fi

# 1D Shockwave files
result=0
filenotfound=0
temp=0
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

if [ $filenotfound -ne 0 ]
then
	if [ $result -ne $(($filenotfound * 2)) ]
	then
		echo "$filenotfound 1DShock Files not found"
	else
		echo "1DShock Files not found"
	fi
else
	if [ $result -eq 0 ]
	then
		echo "1DShock Files are the same"
	else
		echo "1DShock Files are different"
	fi
fi

# 1D Rarefraction files
result=0
filenotfound=0
temp=0
for file in *1d*till5*
do
	diff -q $file ./DiffableResults/1DRare/$file &>/dev/null
	temp=$?
	if [ $temp -eq 2 ]
	then
		filenotfound=$(($filenotfound + 1))
	fi
	result=$(($result + $temp))
done

if [ $filenotfound -ne 0 ]
then
	if [ $result -ne $(($filenotfound * 2)) ]
	then
		echo "$filenotfound 1DRare Files not found"	
	else
		echo "1DRare Files not found"
	fi
else
	if [ $result -eq 0 ]
	then
		echo "1DRare Files are the same"
	else
		echo "1DRare Files are different"
	fi
fi
