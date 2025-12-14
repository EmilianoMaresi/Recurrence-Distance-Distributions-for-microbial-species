#!/usr/bin/env bash
set -e

IMAGE_NAME="rdd_tool"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Installing RDD Docker image..."

# ---- Check Docker availability ----
if ! command -v docker >/dev/null 2>&1; then
    echo "ERROR: Docker is not installed or not in PATH."
    exit 1
fi

if ! docker info >/dev/null 2>&1; then
    echo "ERROR: Docker is not running or permission denied."
    exit 1
fi

# ---- Check if image already exists ----
if docker image inspect "$IMAGE_NAME" >/dev/null 2>&1; then
    read -p "Docker image '$IMAGE_NAME' already exists. Reinstall (overwrite)? [y/n]: " REPLY
    case "$REPLY" in
        [Yy]|[Yy][Ee][Ss])
            echo "Rebuilding Docker image..."
            ;;
        *)
            echo "Installation cancelled. Using existing image."
            echo "You can run the tool using:"
            echo "  ./launch_rdd.sh"
            exit 0
            ;;
    esac
fi

# ---- Build image ----
docker build -t "$IMAGE_NAME" "$SCRIPT_DIR"

echo "Installation completed successfully."
echo "You can now run the tool using:"
echo "  ./launch_rdd.sh"

