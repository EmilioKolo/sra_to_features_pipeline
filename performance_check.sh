

# SRR13051242: 19.55 Gbytes, 55.14 Gpb
# SRR13051091: 29.54 Gbytes, 93.59 Gpb
# SRR13051235: 55.69 Gbytes, 167.36 Gpb

SRA_IDS=("SRR13051242" "SRR13051091" "SRR13051235")
WAIT_TIME=5

# Create the output directory
mkdir -p ./metrics_results

# Loop through each SRA ID
for SRA_ID in "${SRA_IDS[@]}"; do
    echo "Starting metrics measuring for SRA ID: $SRA_ID"
    # Define file paths based on the current SRA ID
    METRICS_FILE="./metrics_results/metrics_"$SRA_ID".csv"
    LOG_FILE="./metrics_results/metrics_log_"$SRA_ID".log"

    # Run the pipeline in the background
    sra-pipeline run --sra-id "$SRA_ID" --output-dir ./metrics_results --config config.ini --threads 4 --log-file "$LOG_FILE" &
    PIPELINE_PID=$!

    # Check if the process started successfully
    if ! ps -p "$PIPELINE_PID" > /dev/null; then
        echo "Error: Failed to start pipeline for $SRA_ID. Exiting loop."
        continue
    fi

    # Start monitoring in a loop
    echo "timestamp,cpu_percent,mem_percent,vsz,rss" > "$METRICS_FILE"
    while ps -p "$PIPELINE_PID" > /dev/null; do
        PS_OUTPUT=$(ps -o %cpu,%mem,vsz,rss -p "$PIPELINE_PID" | tail -n 1)
        CURRENT_TIME=$(date +%s)
        echo "$CURRENT_TIME,$PS_OUTPUT" >> "$METRICS_FILE"
        sleep "$WAIT_TIME"
    done

    echo "Pipeline for $SRA_ID completed. Metrics saved to $METRICS_FILE"
    echo "----------------------------------------------------"
done

echo "All pipelines finished."