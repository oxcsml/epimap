#!/bin/bash

set -e

if [ $# -ne 2 ]; then
  echo Usage: monitor-run results_directory job_id
  exit 1
fi

results_directory=$1
jobid=$2
echo Monitoring job "$results_directory"

tmux new-window -n "$results_directory $jobid"
tmux split-window -h
tmux split-window -h
tmux split-window -h
tmux split-window -h
tmux select-pane -t 1
tmux split-window -h
tmux split-window -h
tmux split-window -h
tmux select-pane -t 1
tmux split-window -h
tmux split-window -h
tmux split-window -h
tmux select-layout tiled


for i in {1..10}; do
  tmux select-pane -t $i && tmux send-keys \
    "while true; do \
       tail -n 1000 -f ${results_directory}/output/run_${jobid}_$i.out; \
       sleep 10; \
     done" C-m
done
tmux select-pane -t 11 && tmux send-keys "swatch Rmap" C-m

# && tmux kill-window

