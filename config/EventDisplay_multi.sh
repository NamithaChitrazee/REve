#!/bin/bash

# ==============================
# Mu2e Event Display Multi Script
# ==============================
# This script processes multiple events from a text file
# Each line should contain event numbers in format: subrun,run,event


usage() {
  echo "Usage: $0 --eventfile <file> --dataset <dataset> [--geom <geometry>]
  
  Required:
    --eventfile  : Text file with event numbers (format: run,sub,event per line)
    --dataset    : Dataset file to extract events from
    
  Optional:
    --geom       : Geometry version (default: 2025)
    
  Example:
    ./EventDisplay_multi.sh --eventfile events.txt --dataset mcs.mu2e.ensembleMDS2cMix1BBTriggered.MDC2020ba_best_v1_3.art --geom 2025
  " 1>&2
}

# Function: Exit with error.
exit_abnormal() {
  usage
  exit 1
}

# Initialize variables
EVENTFILE=""
DATASET=""
GEOM=2025

# Parse command-line arguments
while getopts ":-:" options; do
  case "${options}" in
    -)
      case "${OPTARG}" in
        eventfile)
          EVENTFILE=${!OPTIND} OPTIND=$(( $OPTIND + 1 ))
          ;;
        dataset)
          DATASET=${!OPTIND} OPTIND=$(( $OPTIND + 1 ))
          ;;
        geom)
          GEOM=${!OPTIND} OPTIND=$(( $OPTIND + 1 ))
          ;;
        *)
          echo "Error: Unknown option --${OPTARG}"
          exit_abnormal
          ;;
      esac;;
    :)
      echo "Error: -${OPTARG} requires an argument."
      exit_abnormal
      ;;
    *)
      exit_abnormal
      ;;
  esac
done

# Validate required arguments
if [ -z "$EVENTFILE" ] || [ -z "$DATASET" ]; then
  echo "Error: Missing required arguments."
  exit_abnormal
fi

# Check if event file exists
if [ ! -f "$EVENTFILE" ]; then
  echo "Error: Event file '$EVENTFILE' not found."
  exit 1
fi

echo "========== Launching Mu2e/EventDisplay Multi ==============="
echo "Processing events from: $EVENTFILE"
echo "Dataset: $DATASET"
echo "Geometry: $GEOM"

# Create directory for output files
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="EventDisplay_output_${TIMESTAMP}"

echo "Creating output directory: $OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

# Counter for processed events
event_count=0
success_count=0

# Read event file and process each event
while IFS=',' read -r RUN SUBRUN EVENT; do
  # Skip empty lines and comments
  [[ -z "$RUN" || "$RUN" =~ ^# ]] && continue
  
  # Trim whitespace
  RUN=$(echo "$RUN" | xargs)
  SUBRUN=$(echo "$SUBRUN" | xargs)
  EVENT=$(echo "$EVENT" | xargs)
  
  event_count=$((event_count + 1))
  
  echo ""
  echo "--- Processing event $event_count: subrun=$SUBRUN, run=$RUN, event=$EVENT ---"
  
  # Construct the EVENT identifier string in the format "RUN/SUBRUN/EVENT"
  RUN_SUBRUN_EVENT="${RUN}/${SUBRUN}/${EVENT}"
  
  first_three_chars="${DATASET:0:3}"
  
  DATASET_TEMP="${DATASET}"
  
  if [ "$first_three_chars" == "nts" ]; then
    FILE=$(samweb file-lineage parents ${DATASET_TEMP})
    echo "File: ${FILE}"
    NEWNAME=${FILE%.*}
    NEWNAME=${NEWNAME%.*}
    DATASET_TEMP=${NEWNAME}.art
    echo "Dataset name: ${DATASET_TEMP}"
  fi
  
  # Execute the pickEvent command
  echo "Extracting event..."
  pickEvent -e -v "${DATASET_TEMP}" "${RUN_SUBRUN_EVENT}"
  if [ $? -ne 0 ]; then
    echo "Warning: pickEvent command failed for event $event_count (subrun=$SUBRUN, run=$RUN, event=$EVENT)"
    continue
  fi
  
  success_count=$((success_count + 1))
  
done < "$EVENTFILE"

# Move all generated .art files to output directory
echo ""
echo "--- Moving all .art files to $OUTPUT_DIR ---"

ART_FILES_FOUND=0
for art_file in *.art; do
  if [ -f "$art_file" ]; then
    echo "Moving: $art_file"
    mv "$art_file" "$OUTPUT_DIR/"
    ART_FILES_FOUND=$((ART_FILES_FOUND + 1))
  fi
done

if [ $ART_FILES_FOUND -eq 0 ]; then
  echo "Error: No .art files were generated."
  exit 1
fi

# Create a .txt file listing all .art files
echo ""
echo "--- Creating file list ---"

#cd "$OUTPUT_DIR"
FILE_LIST="art_files.txt"

# Generate list of all .art files
ls -1 "$OUTPUT_DIR"/*.art > "$FILE_LIST"
if [ $? -ne 0 ]; then
    echo "Error: Failed to create file list."
    exit 1
fi

echo "Created file list: $FILE_LIST"
cat "$FILE_LIST"

# Run mu2e with the file list and the specified geometry
echo ""
echo "--- Running display on files ---"


mu2e -c EventDisplay/examples/nominal_MDC${GEOM}.fcl -S $FILE_LIST
if [ $? -ne 0 ]; then
    echo "Error: mu2e command failed."
    exit 1
fi

echo "Display finished."

echo ""
echo "========== Processing Complete ==============="
echo "Total events processed: $event_count"
echo "Successful extractions: $success_count"
echo "Art files created: $ART_FILES_FOUND"
echo "Output directory: $OUTPUT_DIR"
echo "File list: $FILE_LIST"

cd -

exit 0
