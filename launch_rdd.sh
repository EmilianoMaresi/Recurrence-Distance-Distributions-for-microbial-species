#!/usr/bin/env bash
set -e

IMAGE_NAME="rdd_tool"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INSTALL_SCRIPT="$SCRIPT_DIR/install.sh"

show_help() {
    cat << EOF
Usage:
  ./launch_rdd.sh download [arguments...]
  ./launch_rdd.sh run [arguments...]

Commands:
  download    Download datasets from NCBI
  run         Run main analysis

Notes:
  - Docker must be installed and running
  - If the Docker image is missing, the script will suggest running ./install.sh

Examples:
  ./example.sh

EOF
}

# ---- Check Docker availability ----
if ! command -v docker >/dev/null 2>&1; then
    echo "ERROR: Docker is not installed or not in PATH."
    exit 1
fi

if ! docker info >/dev/null 2>&1; then
    echo "ERROR: Docker is not running or permission denied."
    exit 1
fi

# ---- Check if image exists, if not prompt to install ----
if ! docker image inspect "$IMAGE_NAME" >/dev/null 2>&1; then
    echo "Docker image '$IMAGE_NAME' not found."
    if [ -f "$INSTALL_SCRIPT" ]; then
        read -p "Do you want to build it now using ./install.sh? [y/n]: " REPLY
        case "$REPLY" in
            [Yy]|[Yy][Ee][Ss])
                bash "$INSTALL_SCRIPT"
                ;;
            *)
                echo "Please run ./install.sh manually before using the tool."
                exit 1
                ;;
        esac
    else
        echo "Install script '$INSTALL_SCRIPT' not found. Cannot proceed."
        exit 1
    fi
fi

# ---- No command provided ----
if [ $# -lt 1 ]; then
    show_help
    exit 1
fi

COMMAND="$1"
shift

case "$COMMAND" in
    download)
        docker run --rm \
            --user "$(id -u):$(id -g)" \
            -v "$SCRIPT_DIR/datasets_lists:/app/datasets_lists" \
            -v "$SCRIPT_DIR/data:/app/data" \
            "$IMAGE_NAME" download_ncbi_dataset.py "$@"
        ;;
    run)
        docker run --rm \
            --user "$(id -u):$(id -g)" \
            -v "$SCRIPT_DIR/data:/app/data" \
            -v "$SCRIPT_DIR/results:/app/results" \
            "$IMAGE_NAME" main.py "$@"
        ;;
    help|-h|--help)
        show_help
        ;;
    *)
        echo "ERROR: Unknown command '$COMMAND'"
        echo
        show_help
        exit 1
        ;;
esac

