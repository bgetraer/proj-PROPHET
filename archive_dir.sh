#!/bin/bash
set -euo pipefail

archive_dir() {
	local DIR="$1"
	local TAR="$DIR/archive.tar"

	# Create the tar if it doesn't exist
	[ -f "$TAR" ] || tar -cf "$TAR" --files-from /dev/null

	# Make a temporary staging folder
	STAGE="$DIR/.staging_archive"
	mkdir -p "$STAGE"


	# Move all files except the tar into staging
	shopt -s dotglob  # include hidden files
	for f in "$DIR"/* "$DIR"/.*; do
		# skip current dir, parent dir, and the tar itself
		[ "$f" = "$DIR/." ] && continue
		[ "$f" = "$DIR/.." ] && continue
		[ "$f" = "$TAR" ] && continue
		[ "$f" = "$STAGE" ] && continue

		echo "checkpoint"
		mv "$f" "$STAGE"
	done
	shopt -u dotglob

	echo "files in temp dir $STAGE"

	# Only proceed if there are files to archive
	if compgen -G "$STAGE/*" > /dev/null; then
		echo "Archiving $DIR → $TAR"
		tar -rf "$TAR" -C "$STAGE" .
		if [ $? -eq 0 ]; then
			echo "Tar successful! Staged files being deleted."
			# Only delete staged files if tar succeeded
			rm -rf "$STAGE"
			echo "Done."
		else
			echo "Tar failed! Staged files not deleted."
		fi
	else
		rmdir "$STAGE" 2>/dev/null
		echo "No new files to archive in $DIR"
	fi
}

archive_dir "$1"
