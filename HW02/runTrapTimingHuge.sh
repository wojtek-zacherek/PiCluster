!#/bin/bash
echo Starting Testing
i=1024
while [ $i -le 1025 ]
do
	mpirun TTRT.bin $i
	i=$(( $i * 2 ))
	echo $i
done
echo Finished Testing
