#!/bin/bash

# Exit on error and undefined variables
set -euo pipefail

# Function to log messages with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Confirm before proceeding
read -p "This will remove all installed components. Are you sure? [y/N] " confirm
if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    log "Cleanup aborted."
    exit 0
fi

# Remove Miniconda installation
if [[ -d "$HOME/miniconda" ]]; then
    log "Removing Miniconda installation..."
    rm -rf "$HOME/miniconda"
    # Also remove conda initialization from shell config files
    for rcfile in ~/.bashrc ~/.bash_profile ~/.zshrc ~/.profile; do
        if [[ -f "$rcfile" ]]; then
            sed -i '/miniconda/d' "$rcfile"
            sed -i '/conda initialize/d' "$rcfile"
        fi
    done
else
    log "Miniconda not found, skipping removal."
fi

# Remove conda environments (if miniconda wasn't removed)
if [[ -d "$HOME/.conda" ]]; then
    log "Removing conda environments and configuration..."
    rm -rf "$HOME/.conda"
fi

# Remove installation logs
if [[ -d "$HOME/install_logs" ]]; then
    log "Removing installation logs..."
    rm -rf "$HOME/install_logs"
fi

# Remove content directory
if [[ -d "$HOME/content" ]]; then
    log "Removing content directory..."
    rm -rf "$HOME/content"
fi

# Remove snpEff installation
if [[ -d "$HOME/content/data/bin/snpEff" ]]; then
    log "Removing snpEff installation..."
    rm -rf "$HOME/content/data/bin/snpEff"
fi

# Remove Kraken2 database
if [[ -d "$HOME/kraken2-db" ]]; then
    log "Removing Kraken2 database..."
    rm -rf "$HOME/kraken2-db"
fi

# Remove Miniconda installer if it exists
if [[ -f "$HOME/miniconda.sh" ]]; then
    log "Removing Miniconda installer..."
    rm -f "$HOME/miniconda.sh"
fi

# Remove CNVpytor data
if [[ -d "$HOME/.cnvpytor" ]]; then
    log "Removing CNVpytor data..."
    rm -rf "$HOME/.cnvpytor"
fi

# Remove Python cache and build files
log "Cleaning Python cache files..."
find "$HOME" -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
find "$HOME" -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true

log "Cleanup completed successfully!"
