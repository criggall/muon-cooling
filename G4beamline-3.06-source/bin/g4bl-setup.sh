#	g4bl-setup.sh - setup PATH to use this instance of G4beamline
#
#	USAGE: source ...path.../g4bl-setup.sh
#	...path... can be absolute or relative.

G4BL_DIR="$(dirname "${BASH_SOURCE[0]}")"
G4BL_DIR="$(\cd "$G4BL_DIR" >&- 2>&- && pwd)"
G4BL_DIR="$(dirname "$G4BL_DIR")"

export PATH="$G4BL_DIR/bin:$PATH"

echo "G4beamline is in '$G4BL_DIR'"

unset G4BL_DIR 
