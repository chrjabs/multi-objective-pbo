{
  description = "Developement environment for MoaMopb.";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs";
    systems.url = "github:nix-systems/default";
  };

  outputs = {
    self,
    nixpkgs,
    systems,
    ...
  }: let
    lib = nixpkgs.lib;
    pkgsFor = lib.genAttrs (import systems) (system: (import nixpkgs {
      inherit system;
      config.allowUnfree = true;
    }));
    forEachSystem = f: lib.genAttrs (import systems) (system: f pkgsFor.${system});
  in {
    devShells = forEachSystem (pkgs: {
      default = let
        cplex =
          pkgs.cplex.override
          {releasePath = /persist/home/christoph/cplex_studio2211.linux_x86_64.bin;};
      in
        pkgs.mkShell rec {
          nativeBuildInputs = with pkgs; [
            julia
            gurobi
            cplex
          ];
          GRB_LICENSE_FILE = "/persist/home/christoph/gurobi.lic";
          CPLEX_STUDIO_BINARIES = "${cplex}/cplex/bin/x86-64_linux";
          LD_LIBRARY_PATH = "${pkgs.gfortran.cc.lib}/lib";
        };
    });
    formatter = forEachSystem (pkgs: pkgs.alejandra);
  };
}
