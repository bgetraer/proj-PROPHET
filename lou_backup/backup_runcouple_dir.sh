#!/usr/bin/env bash
# sync_all.sh
# Run all rsync for ISSM experiments
# Requires: sync_dir.sh in the same directory or on PATH

set -euo pipefail

# --- root directories ---
ROOT_DIR="/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments"
DEST_ROOT_DIR="/u/bgetraer/backup/proj-PROPHET/experiments"

# --- path to sync_dir.sh (adjust if needed) ---
SYNC_SCRIPT="$(dirname "$0")/sync_dir.sh"

# --- SCENARIOS (optional section) ---
declare -a SCENARIO_SOURCE_DIRS=(
"Paris2C/runcoupledRUN02_dt100_ct1296000"
"RCP85/runcoupledRUN02_dt100_ct1296000"
)
declare -a SCENARIO_DEST_DIRS=(
"Paris2C/RUN02/runcoupled"
"RCP85/RUN02/runcoupled"
)

echo "=== Syncing scenario experiments ==="
for i in "${!SCENARIO_SOURCE_DIRS[@]}"; do
	SRC="$ROOT_DIR/${SCENARIO_SOURCE_DIRS[$i]}"
	DST="$DEST_ROOT_DIR/${SCENARIO_DEST_DIRS[$i]}"
	echo "Syncing: $SRC -> lou:$DST"
	"$SYNC_SCRIPT" "$SRC" "$DST"
done

# --- SENSITIVITY EXPERIMENTS ---
echo "=== Syncing sensitivity experiments ==="
SENS_ROOT="$ROOT_DIR/sensitivity_experiments"

# find all directories matching se*/runcoupled
mapfile -t SENS_DIRS < <(find "$SENS_ROOT" -mindepth 2 -maxdepth 2 -type d -name "runcoupled" | sort)

for SRC in "${SENS_DIRS[@]}"; do
	# compute relative path after ROOT_DIR
	REL_PATH="${SRC#$ROOT_DIR/}"
	DST="$DEST_ROOT_DIR/$REL_PATH"
	echo "Syncing: $SRC -> lou:$DST"
	"$SYNC_SCRIPT" "$SRC" "$DST"
done

echo "=== All syncs complete ==="

