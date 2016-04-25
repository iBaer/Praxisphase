#! /bin/bash

result=0
for file in *5till5*
do
	diff -q $file ./DiffableResults/$file
	result=$(($result + $?))
done


if [ $result -eq 0 ]
then
        echo "Files are the same"
else
        echo "files are different"
        echo "$result"
fi


# Temporary check if 1D warefraction wave and one shockwave results are different

result=0
for file in ./DiffableResults/1DRare/*till5*
do
	diff -q $file ./DiffableResults/1DShock/$(basename $file)
	result=$(($result + $?))
done


if [ $result -eq 0 ]
then
        echo "Files are the same"
else
        echo "files are different"
        echo "$result"
fi
