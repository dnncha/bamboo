//! Generate committed test fixtures under `tests/data/`.

use bamboo_noodles::fixtures::{write_tiny_bam, write_tiny_bam_index};
use std::env;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("tests/data"));

    std::fs::create_dir_all(&root)?;

    let bam_path = root.join("tiny.bam");
    write_tiny_bam(&bam_path)?;
    write_tiny_bam_index(&bam_path)?;

    println!("Wrote {} and {}", bam_path.display(), bam_path.with_extension("bam.bai").display());
    Ok(())
}