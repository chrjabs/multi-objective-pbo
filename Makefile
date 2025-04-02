.PHONY: all scuttle moco-openwbo native-pbo branch-and-bound moa-mopb bin-dir fgt-converter gbmosplit

all: scuttle moco-openwbo native-pbo branch-and-bound moa-mopb fgt-converter gbmosplit

scuttle: bin-dir
	cargo +nightly build --release --manifest-path source-code/scuttle/Cargo.toml --target-dir source-code/scuttle/target
	cp source-code/scuttle/target/release/scuttle bin/

moco-openwbo: bin-dir
	cd source-code/moco-openwbo && make
	cp source-code/moco-openwbo/open-wbo bin/

native-pbo: bin-dir
	cmake -S source-code/native-pbo -B source-code/native-pbo/build -DCMAKE_BUILD_TYPE=Release
	cmake --build source-code/native-pbo/build
	cp source-code/native-pbo/build/roundingsat bin/

branch-and-bound: bin-dir
	cmake -S source-code/branch-and-bound -B source-code/branch-and-bound/build -DCMAKE_BUILD_TYPE=Release
	cmake --build source-code/branch-and-bound/build
	cp source-code/branch-and-bound/build/forget21 bin/

moa-mopb:
	cp source-code/MoaMopb/tiny-inst.mopb .
	julia --project=source-code/MoaMopb -e 'using PackageCompiler; create_sysimage(["MoaMopb"], sysimage_path="bin/MoaMopb.so", precompile_execution_file="source-code/MoaMopb/precompile.jl")'
	rm tiny-inst.mopb

fgt-converter: bin-dir
	cargo build --release --manifest-path source-code/scuttle/rustsat/tools/Cargo.toml --target-dir source-code/scuttle/target --bin=mo2ilp --no-default-features --features=cadical
	cp source-code/scuttle/target/release/mo2ilp bin/

gbmosplit: bin-dir
	cargo build --release --manifest-path source-code/scuttle/rustsat/tools/Cargo.toml --target-dir source-code/scuttle/target --bin=gbmosplit --no-default-features --features=cadical
	cp source-code/scuttle/target/release/gbmosplit bin/

bin-dir:
	mkdir -p bin/
