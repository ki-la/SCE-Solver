fn main() -> miette::Result<()> {
    let path = std::path::PathBuf::from("src"); // include path
    let path2 = std::path::PathBuf::from("lib");
    let mut b = autocxx_build::Builder::new("src/min_cut.rs", &[&path, &path2])
        .extra_clang_args(&["-std=c++2a"])
        .build()?;


    b.flag_if_supported("-std=c++17")
        .flag_if_supported("-O3")
        .flag_if_supported("-DNDEBUG")
        .compile("autocxx-min_cut");
    println!("cargo:rerun-if-changed=src/min_cut.rs");

    Ok(())
}
