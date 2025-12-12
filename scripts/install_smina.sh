#!/bin/bash

# Install Smina for Consensus Docking
# Smina is a fork of AutoDock Vina with improved scoring

echo "Installing Smina for Consensus Docking..."
echo ""

# Check system
OS=$(uname -s)
ARCH=$(uname -m)

echo "System: $OS $ARCH"
echo ""

if [[ "$OS" == "Darwin" ]]; then
    # macOS
    echo "macOS detected"

    if [[ "$ARCH" == "arm64" ]]; then
        # Apple Silicon (M1/M2/M3)
        echo "⚠️  Apple Silicon detected"
        echo ""
        echo "Smina pre-built binaries are not available for Apple Silicon (ARM64)."
        echo ""
        echo "OPTIONS:"
        echo ""
        echo "1. Use Rosetta 2 (emulation) - Easier:"
        echo "   Install Homebrew under Rosetta:"
        echo "   arch -x86_64 /bin/bash -c \"\$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\""
        echo "   Then try: arch -x86_64 brew install smina"
        echo ""
        echo "2. Compile from source - Advanced:"
        echo "   git clone https://github.com/mwojcikowski/smina.git"
        echo "   cd smina"
        echo "   mkdir build && cd build"
        echo "   cmake .."
        echo "   make"
        echo "   sudo cp smina /usr/local/bin/"
        echo ""
        echo "3. Use Docker - Recommended for production:"
        echo "   docker pull ghcr.io/mwojcikowski/smina:latest"
        echo ""
        echo "4. Use consensus docking with Vina only:"
        echo "   Your platform will work with single-method validation until Smina is available."
        echo ""
    else
        # Intel Mac
        echo "Intel Mac detected"
        echo "Attempting to download pre-built binary..."

        SMINA_URL="https://sourceforge.net/projects/smina/files/smina.osx/download"

        curl -L "$SMINA_URL" -o smina
        chmod +x smina
        sudo mv smina /usr/local/bin/

        echo "✓ Smina installed to /usr/local/bin/smina"
    fi

elif [[ "$OS" == "Linux" ]]; then
    # Linux
    echo "Linux detected"

    if [[ "$ARCH" == "x86_64" ]]; then
        echo "Downloading pre-built binary..."

        SMINA_URL="https://sourceforge.net/projects/smina/files/smina.static/download"

        wget -O smina "$SMINA_URL"
        chmod +x smina
        sudo mv smina /usr/local/bin/

        echo "✓ Smina installed to /usr/local/bin/smina"
    else
        echo "⚠️  Architecture $ARCH not supported by pre-built binaries"
        echo "Please compile from source: https://github.com/mwojcikowski/smina"
    fi

else
    echo "⚠️  Unsupported operating system: $OS"
    echo "Please visit: https://sourceforge.net/projects/smina/"
fi

echo ""
echo "Testing installation..."
if command -v smina &> /dev/null; then
    echo "✓ Smina is available!"
    smina --version
else
    echo "⚠️  Smina not found in PATH"
    echo ""
    echo "FALLBACK: The platform will work with Vina-only docking."
    echo "Consensus docking requires multiple methods but can run with just Vina."
fi

echo ""
echo "Done!"
