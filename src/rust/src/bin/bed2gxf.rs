use bed2gtf;
use clap::Parser;

use rust::read;

#[derive(clap::ValueEnum, Clone, Debug)]
enum Gxf {
    Gtf,
    Gff
}

#[derive(Parser, Debug)]
#[command(version = "1.0", about, long_about = None)]
/// Converts file
struct Args {
    /// Input Bed file
    #[arg(long, short = 'b', default_value_t = String::from("stdin"))]
    input: String,

    /// Output format; should be either 'gff' or 'gtf'
    #[arg(long, short = 'f', value_enum, default_value_t = Gxf::Gtf)]
    format: Gxf,

    /// Isoforms file; a two-column tab-separated table of gene-to-transcript correspondence
    #[arg(long, short = 'i')]
    isoforms: Option<String>,

    /// Output file to save the results to
    #[arg(long, short = 'o', default_value_t = String::from("stdout"))]
    output: String
}

fn main() {
    let args = Args::parse();
    // read the input Bed file
}