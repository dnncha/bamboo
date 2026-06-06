//! Generate committed test fixtures under `tests/data/`.

use bamboo_noodles::fixtures::{
    write_tiny_bam, write_tiny_bam_index, write_tiny_vcf, write_tiny_vcf_gz, write_tiny_vcf_index,
};
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

    let vcf_path = root.join("tiny.vcf");
    write_tiny_vcf(&vcf_path)?;

    let vcf_gz_path = root.join("tiny.vcf.gz");
    write_tiny_vcf_gz(&vcf_gz_path)?;
    write_tiny_vcf_index(&vcf_gz_path)?;

    println!(
        "Wrote {}, {}, {}, {}, and {}",
        bam_path.display(),
        bam_path.with_extension("bam.bai").display(),
        vcf_path.display(),
        vcf_gz_path.display(),
        format!("{}.tbi", vcf_gz_path.display())
    );
    Ok(())
}