logdir=$1

for log in $logdir/*
do 
    if [ ! -s $log ]; then
        rm $log;
    fi
done