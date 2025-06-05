{
  description = "Python 3.10 environment using uv and pyproject.toml for local project";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.05";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils, ... }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };
        uv = pkgs.uv;

        projectName = "ttdemo";
        venvDir = ".venv";
        pythonVersion = "3.12";

      in
      {
        devShells.default = pkgs.mkShell {
          packages = [ uv ];

          shellHook = ''
            export PROJECT_NAME="${projectName}"
            export VIRTUAL_ENV_DIR="${venvDir}"
            export PYTHON_VERSION="${pythonVersion}"

            # Ensure the specified Python version is installed by uv
            uv python install "$PYTHON_VERSION" || exit 1

            if [ ! -f "pyproject.toml" ]; then
              uv init "$PROJECT_NAME" || exit 1
              echo "A new pyproject.toml was created. Remember to configure [build-system] and [tool.uv.workspace]."
            fi

            if [ ! -d "$VIRTUAL_ENV_DIR" ]; then
              uv venv "$VIRTUAL_ENV_DIR" --python="$PYTHON_VERSION" || exit 1
            fi
            source "$VIRTUAL_ENV_DIR"/bin/activate || exit 1

            uv sync || exit 1
          '';
        };
      }
    );
}
