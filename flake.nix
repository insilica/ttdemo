{
  description = "Python 3.10 environment for ttdemo with rdkit and related packages";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.05";

    outputs = { self, nixpkgs, ... }: {
    devShells."x86_64-linux".default = nixpkgs.legacyPackages."x86_64-linux".mkShell {
        buildInputs = with nixpkgs.legacyPackages."x86_64-linux".python310Packages; [
        rdkit
        pandas
        rapidfuzz
        aiohttp
        tenacity
        pyarrow
        fastparquet
        seaborn
        matplotlib
        scikit-learn
        cython
        blosc2
        ];

        shellHook = ''
        echo "Installing pip-only packages..."
        python3 -m pip install --quiet --no-cache-dir pubchempy FuzzyTM fastcluster
        '';
        };
    };

}
