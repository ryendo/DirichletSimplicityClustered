#!/bin/bash

n=$1

echo matlab -nodisplay -nosplash -nodesktop -r "try; my_intlab_mode_config(${n}); ${FUNC_SCRIPT}(${j}); catch ME; disp(getReport(ME, 'extended')); exit(1); end; exit(0);" 
echo $0
echo $BASH_VERSION