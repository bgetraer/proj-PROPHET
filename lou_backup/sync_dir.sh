#!/usr/bin/env bash
# sync_dir.sh
# Usage: ./sync_dir.sh /path/to/source_dir /remote/destination_dir

set -euo pipefail

SOURCE_DIR="$1"
DEST_DIR="$2"
SYNC_LIST="/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET/lou_backup/sync_list.txt"

mkdir -p "$(dirname "$SYNC_LIST")"
> "$SYNC_LIST"  # clear sync list

# --- Step 1: Collect timestep numbers from issmDiag*.mat ---
mapfile -t T1 < <(find "$SOURCE_DIR" -maxdepth 1 -type f -name "issmDiag*.mat" \
	| sed -E 's/.*issmDiag\.([0-9]+)\.mat/\1/' | sort -n | uniq)

# --- Step 2: Collect timestep numbers from *.save.*.* ---
mapfile -t T2 < <(find "$SOURCE_DIR" -maxdepth 1 -type f -name "*.save.*.*" \
	| sed -E 's/.*\.save\.([0-9]+)\..*/\1/' | sort -n | uniq)

# --- Step 3: Combine and deduplicate ---
mapfile -t T < <(printf "%s\n" "${T1[@]}" "${T2[@]}" | sort -n | uniq)
echo "${#T[@]} timesteps found to sync"

# --- Step 4: Build sync list for all but last timestep ---
N=$(( ${#T[@]} - 1 ))
NFILES=0

cd "$SOURCE_DIR"  # ensure relative paths
for ((i=0; i<N; i++)); do
    t="${T[$i]}"
    # force t to be treated as a decimal string, not a number
    pad_t=$(printf "%010s" "$t" | tr ' ' '0')
    matches=$(find . -maxdepth 1 -type f -name "*.${pad_t}.*" | sed 's|^\./||')
    if [[ -n "$matches" ]]; then
        echo "$matches" >> "$SYNC_LIST"
        count=$(echo "$matches" | wc -l)
        ((NFILES+=count))
    fi
done

echo "$NFILES files found to sync"

# --- Step 5: Dry-run sync ---
echo "Syncing..."
rsync -av --remove-source-files --files-from="$SYNC_LIST" "$SOURCE_DIR"/ "lou:$DEST_DIR"

