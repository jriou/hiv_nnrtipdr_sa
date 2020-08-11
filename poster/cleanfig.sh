#!/bin/bash

for image_file in $(ls figures/)
do
if grep $image_file *.log -c > 1
then
	echo "File $image_file is in use."
else
	echo "File $image_file is not in use."
	mv "figures/$image_file" "figures/moved.$image_file"
fi
done
