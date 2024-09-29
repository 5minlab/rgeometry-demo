## prerequisites

- [rust](https://www.rust-lang.org/tools/install)
- [wasm-pack](https://rustwasm.github.io/wasm-pack/installer/)

```sh
cargo install wasm-pack --version 0.12.1
```

## gui debugger

```sh
# run on a gui environment
cargo run --release -p gui
```

## build wasm artifacts
```sh
export TARGET_DIR= # target directory
(cd wasm && wasm-pack build -t web -d "${TARGET_DIR}/vendor/asm")
(cd wasm && wasm-pack build -t nodejs -d "${TARGET_DIR}/vendor/asm-node")
```
