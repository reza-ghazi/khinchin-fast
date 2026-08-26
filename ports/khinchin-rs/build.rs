fn main() {
    cc::Build::new().file("src/shim.c").compile("khinchin_shim");
    println!("cargo:rustc-link-lib=flint");
    println!("cargo:rerun-if-changed=src/shim.c");
}
