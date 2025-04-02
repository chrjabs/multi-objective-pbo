{
  description = "RoundingSAT";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs";
    nixpkgs.follows = "nix-config/nixpkgs";
    systems.url = "github:nix-systems/default-linux";
    nix-config.url = "github:chrjabs/nix-config";
  };

  outputs = {
    self,
    nixpkgs,
    systems,
    nix-config,
  }: let
    lib = nixpkgs.lib;
    pkgsFor = lib.genAttrs (import systems) (system: (import nixpkgs {
      inherit system;
      overlays = builtins.attrValues nix-config.overlays;
    }));
    forEachSystem = f: lib.genAttrs (import systems) (system: f pkgsFor.${system});
  in {
    devShells = forEachSystem (pkgs: {
      default = pkgs.mkShell.override {stdenv = pkgs.clangStdenv;} {
        nativeBuildInputs = with pkgs; [
          cmake
          boost
          lldb
          veripb
        ];
      };
    });
  };
}
