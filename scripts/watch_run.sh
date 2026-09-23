#!/bin/bash
# Keep a detached solver run alive until it reaches `done`.
#
# Why this exists: the run is hours long and this container is reclaimed on its
# own schedule.  The solver checkpoints every 180 s, so a death costs at most
# three minutes -- but only if something relaunches it.  Three runs were lost
# that way earlier in this project, at 22%, 22% and 2%.
#
# It relaunches the SAME command, which resumes from the checkpoint: the
# settings tag inside the checkpoint makes a resume under different physics
# impossible (R-2), so a relaunch is safe by construction rather than by care.
#
#   scripts/watch_run.sh <status.json> <log> <command-file>
set -u
STATUS="$1"; LOG="$2"; CMDFILE="$3"
JOURNAL="${STATUS%.status.json}.watch.log"

note() { echo "[$(date -u +%FT%TZ)] $*" >> "$JOURNAL"; }

note "watchdog started"
restarts=0
while true; do
    state=$(python3 -c "
import json,sys
try: print(json.load(open('$STATUS'))['state'])
except Exception: print('unknown')
" 2>/dev/null)

    case "$state" in
        done)   note "run finished cleanly"; exit 0 ;;
        failed) note "run reported failure; not relaunching"; exit 1 ;;
    esac

    pid=$(python3 -c "
import json
try: print(json.load(open('$STATUS')).get('pid',0))
except Exception: print(0)
" 2>/dev/null)

    if [ "${pid:-0}" -gt 0 ] && kill -0 "$pid" 2>/dev/null; then
        sleep 60
        continue
    fi

    # The process is gone but the status says it was still working.  Resume.
    restarts=$((restarts + 1))
    if [ "$restarts" -gt 20 ]; then
        note "giving up after $restarts relaunches"; exit 1
    fi
    note "process $pid is gone at state=$state; relaunch #$restarts"
    cd "$(dirname "$0")/.." || exit 1
    setsid nohup env PYTHONPATH=. bash "$CMDFILE" >> "$LOG" 2>&1 &
    disown
    sleep 90
done
