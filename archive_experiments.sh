#!/usr/bin/env bash
set -euo pipefail

archive_experiments(){
	local LOG_DIR="/u/bgetraer/backup/proj-PROPHET/experiments/archive_logs"
	mkdir -p "$LOG_DIR"
	local LOG_FILE="$LOG_DIR/archive_$(date +%Y%m%d_%H%M%S).log"

	local ALL_DIRS=(
	"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se002_KN_constant_clim/runcoupled"
	"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se003_KN_monthly_clim/runcoupled"
	"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se004_PW700_constant/runcoupled"
	"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se005_PW600_constant/runcoupled"
	"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se006_PW800_constant/runcoupled"
	"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se007_PW700_amp50_per2/runcoupled"
	"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se008_PW700_amp50_per5/runcoupled"
	"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se009_PW700_amp50_per10/runcoupled"
	"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se010_PW700_amp100_per2/runcoupled"
	"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se011_PW700_amp100_per5/runcoupled"
	"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se012_PW700_amp100_per10/runcoupled"
	)
	#local ALL_DIRS=(
	#"/u/bgetraer/backup/proj-PROPHET/experiments/Paris2C/RUN02/runcoupled"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/RCP85/RUN02/runcoupled"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se002_KN_constant_clim/runcoupled"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se003_KN_monthly_clim/runcoupled"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se004_PW700_constant/runcoupled"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se005_PW600_constant/runcoupled"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se006_PW800_constant/runcoupled"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se007_PW700_amp50_per2/runcoupled"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se008_PW700_amp50_per5/runcoupled"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se009_PW700_amp50_per10/runcoupled"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se010_PW700_amp100_per2/runcoupled"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se011_PW700_amp100_per5/runcoupled"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se012_PW700_amp100_per10/runcoupled"
	#)
	#local ALL_DIRS=(
	#"/u/bgetraer/backup/proj-PROPHET/experiments/Paris2C/RUN02/runocean"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/RCP85/RUN02/runocean"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se002_KN_constant_clim/runocean"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se003_KN_monthly_clim/runocean"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se004_PW700_constant/runocean"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se005_PW600_constant/runocean"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se006_PW800_constant/runocean"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se007_PW700_amp50_per2/runocean"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se008_PW700_amp50_per5/runocean"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se009_PW700_amp50_per10/runocean"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se010_PW700_amp100_per2/runocean"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se011_PW700_amp100_per5/runocean"
	#"/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se012_PW700_amp100_per10/runocean"
	#)


	for DIR in "${ALL_DIRS[@]}"; do
		if [ -d "$DIR" ]; then
			echo "Processing: $DIR"
			"/u/bgetraer/backup/proj-PROPHET/scripts/archive_dir.sh" "$DIR" >> "$LOG_FILE" 2>&1
		else
			echo "Directory not found: $DIR" | tee -a "$LOG_FILE"
		fi
	done

	echo "All experiments processed. Log saved to $LOG_FILE"
}

archive_experiments
