#!/bin/bash

TARGET_PID=$1
OUTPUT_FILE="system_memory_usage_new.csv"

if [ -f $OUTPUT_FILE ]; then
    rm -f $OUTPUT_FILE
fi

# Function to get all descendant PIDs of a given PID
get_descendant_pids() {
    local parent_pid=$1
    echo $parent_pid
    local children=$(pgrep -P $parent_pid)
    for child in $children; do
        get_descendant_pids $child
    done
}

# Function to measure CPU memory usage of all PIDs using Pss
measure_cpu_memory_usage() {
    local pids=("$@")
    local total_mem=0
    for pid in "${pids[@]}"; do
        if [ -e /proc/$pid/smaps ]; then
            # Sum up the Pss values from the smaps file
            local mem=$(awk '/Pss/ {sum += $2} END {print sum}' /proc/$pid/smaps)
            total_mem=$((total_mem + mem))
        fi
    done
    total_mem2=$(echo "scale=6; $total_mem / 1024 / 1024" | bc) # Convert to GB
    echo $total_mem2
}


# Function to measure GPU memory usage for the given PIDs
measure_gpu_memory_usage() {
    local pids=("$@")
    local total_gpu_mem=0
    local gpu_mem_info=$(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader,nounits)
    for pid in "${pids[@]}"; do
        local gpu_mem=$(echo "$gpu_mem_info" | awk -F',' -v pid="$pid" '$1 == pid {sum+=$2} END {print sum}')
        total_gpu_mem=$((total_gpu_mem + gpu_mem))
    done
    total_gpu_mem2=$(echo "scale=6; $total_gpu_mem / 1024" | bc)
    echo $total_gpu_mem2
}
# Log the initial memory usage as 0
echo "$(date '+%Y-%m-%d %H:%M:%S.%N') 0 GB 0 GB" >> $OUTPUT_FILE
sleep 1

# Main monitoring loop
while true; do
    start_time=$(date +%s.%N)

    if ! ps -p $TARGET_PID > /dev/null 2>&1; then
        echo "Target process $TARGET_PID has terminated."
        break
    fi

    # Get all descendant PIDs
    all_pids=($(get_descendant_pids $TARGET_PID))

    # Measure total CPU memory usage
    total_cpu_memory=$(measure_cpu_memory_usage "${all_pids[@]}")

    # Measure total GPU memory usage
    total_gpu_memory=$(measure_gpu_memory_usage "${all_pids[@]}")

    # Log the result with a timestamp
    echo "$(date '+%Y-%m-%d %H:%M:%S.%N') $total_cpu_memory GB $total_gpu_memory GB" >> $OUTPUT_FILE

    # Check if all processes have terminated
    all_terminated=true
    for pid in "${all_pids[@]}"; do
        if ps -p $pid > /dev/null 2>&1; then
            all_terminated=false
            break
        fi
    done

    if $all_terminated; then
        echo "Monitoring completed. All processes have terminated."
        break
    fi

    # Calculate elapsed time and sleep for the remainder of 1 second
    end_time=$(date +%s.%N)
    elapsed_time=$(echo "$end_time - $start_time" | bc)
    sleep_time=$(echo "1 - $elapsed_time" | bc)
    if (( $(echo "$sleep_time > 0" | bc -l) )); then
        sleep $sleep_time
    fi
done

echo "Monitoring completed. All processes have terminated."