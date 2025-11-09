#!/bin/bash

# Monitor memory usage of a target PID and all descendants
# Samples exactly every 1s (aligned to a fixed schedule) and logs with consistent formatting.

TARGET_PID=$1
STAGE_FILE=${2:-}
GPU_TOGGLE_RAW=${3:-}
OUTPUT_FILE="system_memory_usage_new.csv"

if [ -z "$TARGET_PID" ]; then
    echo "Usage: $0 <target_pid>" >&2
    exit 1
fi

if [ -f "$OUTPUT_FILE" ]; then
    rm -f "$OUTPUT_FILE"
fi

# Get all descendant PIDs of a given PID (recursive)
# Collect descendant PIDs (including the root) in one pass using ps (faster, fewer races)
collect_descendant_pids() {
    local root_pid=$1
    # If ps fails, just return root
    ps -eo pid=,ppid= 2>/dev/null | awk -v root="$root_pid" '
        {
            # Build children list for each ppid
            children[$2] = (children[$2] " " $1)
        }
        END {
            # BFS from root
            qn = 1; qi = 1; q[qn] = root; seen[root] = 1;
            while (qi <= qn) {
                cur = q[qi++];
                if (cur != "") print cur;
                n = split(children[cur], arr, " ");
                for (i = 1; i <= n; i++) {
                    c = arr[i];
                    if (c != "" && !(c in seen)) { seen[c] = 1; q[++qn] = c; }
                }
            }
        }
    ' || echo "$root_pid"
}

"${SHELL:-/bin/bash}" -c ':' >/dev/null 2>&1 # ensure $SHELL set

# Optional faster scope: use PGID to capture process group members instead of full tree
SCOPE_MODE=${MEMORY_MONITOR_SCOPE:-tree}
SCOPE_MODE=$(echo "$SCOPE_MODE" | tr '[:upper:]' '[:lower:]')
ROOT_PGID=$(ps -o pgid= -p "$TARGET_PID" 2>/dev/null | awk '{print $1}' )
collect_pgid_pids() {
    if [ -n "$ROOT_PGID" ]; then
        ps -o pid= -g "$ROOT_PGID" 2>/dev/null | awk '{print $1}'
    else
        echo "$TARGET_PID"
    fi
}

# Measure total CPU (system RAM) usage
# mode rss (default): fast, sums RSS via one or few ps calls
# mode pss: accurate shared-accounting via smaps_rollup/smaps (heavier)
measure_cpu_memory_usage_kb() {
    local pids=("$@")
    local mode="${MEMORY_MONITOR_MODE:-pss}"
    mode=$(echo "$mode" | tr '[:upper:]' '[:lower:]')
    if [ "$mode" = "rss" ]; then
        local total_kb=0
        local chunk=500
        local start=0
        local end=0
        local n=${#pids[@]}
        while [ $start -lt $n ]; do
            end=$(( start + chunk ))
            if [ $end -gt $n ]; then end=$n; fi
            local sub=("${pids[@]:$start:$((end-start))}")
            # Build comma-separated PID list
            local list=""
            local i
            for i in "${sub[@]}"; do
                [ -n "$i" ] && list+="$i,"
            done
            list=${list%,}
            if [ -n "$list" ]; then
                # Sum RSS in kB across this chunk
                local sum_chunk
                sum_chunk=$(ps -o rss= -p "$list" 2>/dev/null | awk '{s+=int($1)} END{print s+0}')
                total_kb=$(( total_kb + ${sum_chunk:-0} ))
            fi
            start=$end
        done
        echo "$total_kb"
        return
    fi
    # PSS path
    local total_kb=0
    local mem_kb=0
    local pid
    for pid in "${pids[@]}"; do
        if [ -r "/proc/$pid/smaps_rollup" ]; then
            mem_kb=$(awk '/^Pss:/ {print $2; exit}' "/proc/$pid/smaps_rollup" 2>/dev/null || echo 0)
        elif [ -r "/proc/$pid/smaps" ]; then
            mem_kb=$(awk '/^Pss:/ {sum += $2} END {print sum+0}' "/proc/$pid/smaps" 2>/dev/null || echo 0)
        else
            mem_kb=$(ps -o rss= -p "$pid" 2>/dev/null | awk '{print int($1)+0}' || echo 0)
        fi
        total_kb=$(( total_kb + mem_kb ))
    done
    echo "$total_kb"
}

# Resolve GPU toggle once for reuse
GPU_TOGGLE_RESOLVED=${GPU_TOGGLE_RAW}
if [ -z "$GPU_TOGGLE_RESOLVED" ]; then
    GPU_TOGGLE_RESOLVED="${MEMORY_MONITOR_GPU:-off}"
fi
GPU_TOGGLE_RESOLVED=$(echo "$GPU_TOGGLE_RESOLVED" | tr '[:upper:]' '[:lower:]')
GPU_ENABLED=1
if [ "$GPU_TOGGLE_RESOLVED" = "off" ] || [ "$GPU_TOGGLE_RESOLVED" = "0" ] || [ "$GPU_TOGGLE_RESOLVED" = "no" ]; then
    GPU_ENABLED=0
fi

# Measure total GPU memory usage (MiB) for the given PIDs using procfs if available, else fall back to nvidia-smi
measure_gpu_memory_usage_mib() {
    local pids=("$@")
    local total_mib=0
    # Check toggle: skip work entirely if disabled
    if [ $GPU_ENABLED -eq 0 ]; then
        echo 0
        return
    fi

    # Try fast procfs path first
    local have_proc=0
    for f in /proc/driver/nvidia/gpus/*/processes; do
        [ -r "$f" ] || continue
        have_proc=1
        # Sum MiB for our PIDs found in this file
        local add
        add=$(awk -v PIDS="$(printf '%s ' "${pids[@]}")" '
            BEGIN {
                split(PIDS, arr)
                for (i in arr) if (arr[i] != "") want[arr[i]] = 1
                pid = ""; used = 0; sum = 0
            }
            /Process ID:/ { pid=$3 }
            /Used GPU Memory:/ {
                # Handle formats like: "Used GPU Memory: 123 MiB" or "Memory: 123 MiB"
                for (i=1; i<=NF; i++) if ($i ~ /^[0-9]+$/) { used=$i; break }
                unit=$(NF)
                mib = (unit ~ /GiB/) ? used*1024 : (unit ~ /KiB/) ? used/1024 : used
                if (pid in want) sum += mib
                pid = ""; used = 0
            }
            END { print sum+0 }
        ' "$f")
        total_mib=$(( total_mib + ${add:-0} ))
    done

    if [ $have_proc -eq 0 ]; then
        if command -v nvidia-smi >/dev/null 2>&1; then
            local gpu_mem_info pid gpu_mem
            gpu_mem_info=$(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader,nounits 2>/dev/null || true)
            if [ -n "$gpu_mem_info" ]; then
                for pid in "${pids[@]}"; do
                    gpu_mem=$(awk -F',' -v pid="$pid" '$1 == pid {sum+=$2} END {print sum+0}' <<<"$gpu_mem_info")
                    total_mib=$(( total_mib + gpu_mem ))
                done
            fi
        fi
    fi

    echo "$total_mib"
}

# Convert kB -> GiB with 6 decimals
kb_to_gib() {
    awk -v kb="$1" 'BEGIN { printf "%.6f", kb/1048576 }'
}

# Convert MiB -> GiB with 6 decimals
mib_to_gib() {
    awk -v mib="$1" 'BEGIN { printf "%.6f", mib/1024 }'
}

# Align sampling to a fixed 1s schedule
ns_now() { date +%s%N; }

interval_ns=1000000000
now_ns=$(ns_now)
# First tick: next whole second from now
next_ns=$(( (now_ns/interval_ns + 1) * interval_ns ))

# This monitor's own PID (exclude from sums)
SELF_PID=$$
exclude_self() {
    local out=()
    local pid
    for pid in "$@"; do
        if [ -n "$pid" ] && [ "$pid" != "$SELF_PID" ]; then
            out+=("$pid")
        fi
    done
    printf "%s\n" "${out[@]}"
}

# Build PID list = TARGET + (PGID or tree) minus SELF, deduped
gather_pids() {
    local list
    if [ "$SCOPE_MODE" = "pgid" ]; then
        list=$(collect_pgid_pids)
    else
        list=$(collect_descendant_pids "$TARGET_PID")
    fi
    printf "%s\n%s\n" "$TARGET_PID" "$list" | awk -v self="$SELF_PID" 'NF && $0!=self { if(!seen[$0]++) print }'
}

# Main monitoring loop
while true; do
    # Stop if target process has ended and no descendants remain
    if ! ps -p "$TARGET_PID" > /dev/null 2>&1; then
        # double-check descendants; if none, finish
        mapfile -t all_pids < <(gather_pids)
        if [ "${#all_pids[@]}" -eq 0 ]; then
            break
        fi
    fi

    # Sleep until scheduled tick
    now_ns=$(ns_now)
    sleep_ns=$(( next_ns - now_ns ))
    if [ $sleep_ns -gt 0 ]; then
        sleep_sec=$(awk -v ns="$sleep_ns" 'BEGIN { printf "%.9f", ns/1000000000 }')
        sleep "$sleep_sec"
    else
        # If we are late, fast-forward the schedule to the next future tick
        missed=$(( (-sleep_ns)/interval_ns + 1 ))
        next_ns=$(( next_ns + missed*interval_ns ))
    fi

    # Sample exactly once per scheduled second
    mapfile -t all_pids < <(gather_pids)
    mapfile -t all_pids < <(exclude_self "${all_pids[@]}")

    total_cpu_kb=$(measure_cpu_memory_usage_kb "${all_pids[@]}")
    if [ $GPU_ENABLED -eq 1 ]; then
        total_gpu_mib=$(measure_gpu_memory_usage_mib "${all_pids[@]}")
    else
        total_gpu_mib=0
    fi

    cpu_gib=$(kb_to_gib "$total_cpu_kb")
    if [ $GPU_ENABLED -eq 1 ]; then
        gpu_gib=$(mib_to_gib "$total_gpu_mib")
    fi

    # Read current stage label (optional)
    stage_label=""
    if [ -n "$STAGE_FILE" ] && [ -r "$STAGE_FILE" ]; then
        stage_label=$(head -n1 "$STAGE_FILE" 2>/dev/null | tr '\n' ' ')
    fi

    # Log the scheduled tick timestamp (not the completion time) for stable 1s cadence
    tick_ns=$(( next_ns - interval_ns ))
    tick_sec=$(( tick_ns / 1000000000 ))
    tick_nano=$(( tick_ns % 1000000000 ))
    ts_date=$(date -d "@${tick_sec}" '+%Y-%m-%d %H:%M:%S')
    if [ $GPU_ENABLED -eq 1 ]; then
        if [ -n "$stage_label" ]; then
            printf "%s.%09d %s GB %s GB STAGE:%s\n" "$ts_date" "$tick_nano" "$cpu_gib" "$gpu_gib" "$stage_label" >> "$OUTPUT_FILE"
        else
            printf "%s.%09d %s GB %s GB\n" "$ts_date" "$tick_nano" "$cpu_gib" "$gpu_gib" >> "$OUTPUT_FILE"
        fi
    else
        if [ -n "$stage_label" ]; then
            printf "%s.%09d %s GB STAGE:%s\n" "$ts_date" "$tick_nano" "$cpu_gib" "$stage_label" >> "$OUTPUT_FILE"
        else
            printf "%s.%09d %s GB\n" "$ts_date" "$tick_nano" "$cpu_gib" >> "$OUTPUT_FILE"
        fi
    fi

    # Schedule next tick (exact 1s cadence)
    next_ns=$(( next_ns + interval_ns ))

    # If all processes are gone, exit after logging this tick
    all_terminated=true
    for pid in "${all_pids[@]}"; do
        if ps -p "$pid" > /dev/null 2>&1; then
            all_terminated=false
            break
        fi
    done
    $all_terminated && break
done

echo "Monitoring completed. All processes have terminated."
