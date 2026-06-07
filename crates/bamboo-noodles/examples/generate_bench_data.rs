//! Generate synthetic benchmark BAMs under `benchmarks/data/`.

use bamboo_noodles::fixtures::{
    write_bench_bam, write_bench_bam_index, write_bench_cram, write_bench_cram_index,
    write_bench_fasta, write_bench_fasta_index,
};
use std::env;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("benchmarks/data"));

    let record_count: usize = env::args()
        .nth(2)
        .and_then(|value| value.parse().ok())
        .unwrap_or(100_000);

    std::fs::create_dir_all(&root)?;

    let bam_path = root.join(format!("bench_{record_count}.bam"));
    write_bench_bam(&bam_path, record_count)?;
    write_bench_bam_index(&bam_path)?;

    let cram_path = root.join(format!("bench_{record_count}.cram"));
    write_bench_cram(&cram_path, record_count)?;
    write_bench_cram_index(&cram_path)?;

    let fasta_path = root.join(format!("bench_{record_count}.fasta"));
    write_bench_fasta(&fasta_path)?;
    write_bench_fasta_index(&fasta_path)?;

    println!(
        "Wrote {}, {}, {}, {}, {}, and {}",
        bam_path.display(),
        bam_path.with_extension("bam.bai").display(),
        cram_path.display(),
        format!("{}.crai", cram_path.display()),
        fasta_path.display(),
        format!("{}.fai", fasta_path.display())
    );
    Ok(())
}