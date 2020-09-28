#!/bin/bash

set -e

if [ $# -ne 1 ]; then
  echo Usage: tmux-monitor job-id
  exit 1
fi

jobid=$1
echo Monitoring job "$jobid"

tmux new-window -n "job_$jobid"
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
    "while true; do tail -n 1000 -f slurm/output/run_${jobid}_$i.out; sleep 10; done" C-m
done
tmux select-pane -t 11 && tmux send-keys "swatch" C-m

# && tmux kill-window

