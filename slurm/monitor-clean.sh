
#!/bin/bash

trap 'echo monitor-clean: Failed before finishing' ERR

if [ $# -ne 2 ]; then
  echo Usage: monitor-clean clean_directory job_id
  exit 1
fi

clean_directory=$1
jobid=$2
echo Monitoring "$clean_directory"

tmux new-window -n "monitor $clean_directory"
tmux split-window -h

tmux select-pane -t 1 && tmux send-keys \
    "while true; do \
       clear; \
       tmux clear-history; \
       for i in {1..348}; do \
         file=${clean_directory}/output/clean_${jobid}_\${i}.out
         if [ -f \$file ]; then
           echo \${i}: `tail -n 1 \${file}`; \
         fi
       done
       sleep 10; \
     done" C-m
tmux select-pane -t 2 && tmux send-keys "swatch" C-m

# && tmux kill-window

