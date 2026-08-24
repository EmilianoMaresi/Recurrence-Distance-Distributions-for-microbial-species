#!/usr/bin/env bash
set -e

IMAGE_NAME="rdd_tool"

echo "Uninstalling RDD Docker image..."

# ---- Check Docker availability ----
if ! command -v docker >/dev/null 2>&1; then
    echo "ERROR: Docker is not installed or not in PATH."
    exit 1
fi

# ---- Check if the image exists ----
if ! docker image inspect "$IMAGE_NAME" >/dev/null 2>&1; then
    echo "Docker image '$IMAGE_NAME' does not exist or is already uninstalled."
    exit 0
fi

# ---- Confirm uninstallation ----
read -p "Are you sure you want to remove the Docker image '$IMAGE_NAME'? [y/n]: " REPLY
case "$REPLY" in
    [Yy]|[Yy][Ee][Ss])
        echo "Removing Docker image..."
        docker rmi "$IMAGE_NAME"
        echo "Uninstallation completed successfully."
        ;;
    *)
        echo "Uninstallation cancelled."
        exit 0
        ;;
esac
