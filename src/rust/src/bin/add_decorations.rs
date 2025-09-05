use clap::Parser;
use cubiculum::extract::extract::parse_bed;
use cubiculum::structs::structs::{BedEntry, Coordinates};
use fxhash::FxHashMap;
use phf::phf_map;
use std::fs::File;
use std::io::{BufRead, BufReader, Write, stdout};
use std::path::Path;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
/// Given a Bed file and a TOGA2 mutation report file, prepare 
/// 
struct Args {

    /// path to the Bed12 file to decorate
    #[arg(long, short = 'b')]
    bed_file: String,

    ///  path to mutation file
    #[arg(long, short = 'm')]
    mutation_file: String,

    /// path to the output file
    #[arg(long, short = 'o')]
    output: String
}

const MUT2SHAPE: phf::Map<&'static str, &'static str> = phf_map!{
    "FS_DEL" => "Triangle",
    "FS_INS" => "InvTriangle",
    "STOP" => "Star",
    "SSMA" => "Square",
    "SSMD" => "SSMD"
};

const ACCEPTOR: &str = "SSMA";
const DONOR: &str = "SSMD";

const DELETED: &str = "Exon is deleted";
const MISSING: &str = "Exon is missing";

const MASKED: &str = "211,211,211,255";
const NOT_MASKED: &str = "198,0,9,255";
const GLYPH: &str = "glyph";
// static MUT2SHAPE: FxHashMap<&str, String> = FxHashMap::from_iter([
//     ("", String::from("")),
// ]);

#[derive(Clone)]
struct PointMutation {
    chrom: String,
    start: u64,
    end: u64,
    exon: String,
    mut_type: String,
    color: &'static str
}

fn main() {
    let args = Args::parse();
    let mut proj2muts: FxHashMap<String, Vec<PointMutation>> = FxHashMap::default();
    // read the mutation file first
    let mutations = {
        let path = File::open(args.mutation_file).unwrap();
        Box::new(BufReader::new(path)) as Box<dyn BufRead>
    };
    for (i, line_) in mutations.lines().enumerate() {
        if let Ok(line) = line_ {
            let mut_comps: Vec<&str> = line.split('\t').collect::<Vec<&str>>();
            if mut_comps.len() == 0 {continue};
            let proj: &str = mut_comps[0];
            // ignore fragmented projections; due to the current naming scheme, 
            // mutations cannot be attributed to fragments unambiguously
            if proj.contains(',') {continue};
            let mut_type: &str = mut_comps[7];
            // ignore mutations outside of the point mutation list§
            if !MUT2SHAPE.contains_key(mut_type) {continue};
            let mask_reason: &str = mut_comps[11];
            // ignore mutations corresponding to deleted and missing exons 
            if mask_reason == DELETED || mask_reason == MISSING {continue};
            let chrom: &str = mut_comps[4];
            let start: u64 = mut_comps[5]
                .parse::<u64>()
                .expect(&format!(
                    "Improper mutation file formatting at line {}: \"start\" field is not a valie integer", i+1
                ));
            let end: u64 = mut_comps[6]
                .parse::<u64>()
                .expect(&format!(
                    "Improper mutation file formatting at line {}: \"end\" field is not a valie integer", i+1
                ));
            let color: &str = match mut_comps[9] {
                "MASKED" => {MASKED},
                "NOT_MASKED" => {NOT_MASKED},
                _ => {
                    panic!(
                        "Improper mutation file formatting at line {}: \"is_masked\" field contains illegal value", i+1
                    )
                }
            };
            let exon: &str = mut_comps[1];
            let mut_entry: PointMutation = PointMutation { 
                chrom: chrom.to_string(), 
                start: start, 
                end: end, 
                exon: exon.to_string(),
                mut_type: mut_type.to_string(), 
                color: color 
            };
            proj2muts.entry(proj.to_string())
                .and_modify(|x| x.push(mut_entry.clone()))
                .or_insert(vec![mut_entry]);
        }
    }

    // create an output file handle
    let mut output_file = match args.output.as_str() {
        "stdout" => {Box::new(stdout()) as Box<dyn Write>},
        _ => {
            let path = Path::new(&args.output);
            Box::new(File::create(&path).unwrap()) as Box<dyn Write>
        }
    };

    //then, read the Bed file 
    let bed_file = {
        let path = File::open(args.bed_file).unwrap();
        Box::new(BufReader::new(path)) as Box<dyn BufRead>
    };
    for line_ in bed_file.lines() {
        if let Ok(line) = line_ {
            let bed_entry: BedEntry = match parse_bed(line, 12, true) {
                Some(entry) => {entry},
                None => {continue}
            };
            let proj: &str = bed_entry.name().unwrap();
            // skip fragmented projections - TOGA2 currently does not support 
            if proj.contains(',') {continue};
            if let Some(muts) = proj2muts.get(proj) {
                // extract projection values needed for the mutation Bed record
                let strand: bool = bed_entry.strand().unwrap();
                let strand_lit: char = if strand {'+'} else {'-'};
                let tr_chrom: &str = bed_entry.chrom().unwrap();
                let tr_start: &u64 = bed_entry.start().unwrap();
                let tr_end: &u64 = bed_entry.end().unwrap();
                let tr_coords: String = format!("{}:{}-{}:{}", tr_chrom, tr_start, tr_end, proj);
                for mutation in muts {
                    // point mutations have a total length of 1; 
                    // meanwhile, mutation files provided affected codon's coordinate;
                    // therefore, one has to first infer the exact coordinate
                    let first_coord: u64 = match mutation.exon.contains('_') {
                        true => {
                            // mutation spans over multiple exons;
                            // acceptors should be mapped to mutation's end,
                            // the rest go to the first coordinate as always
                            if &mutation.mut_type == ACCEPTOR {
                                if strand {mutation.end} else {mutation.start}
                            } else {
                                if strand {mutation.start} else {mutation.end}
                            }
                        },
                        false => {
                            // mutation occurs within one exon;
                            // unless it is a donor mutation, simply place it 
                            // at the beginning of the codon
                            if &mutation.mut_type == DONOR {
                                if strand {mutation.end} else {mutation.start}
                            } else {
                                if strand {mutation.start} else {mutation.end}
                            }
                        }
                    };
                    let last_coord: u64 = first_coord + 1; // TODO: Chromosome size exceeding safeguard?
                    let shape = MUT2SHAPE.get(&mutation.mut_type).unwrap();
                    let mut_name: String = format!(
                        "{}:{}:{}:{}", first_coord, last_coord, mutation.exon, mutation.mut_type
                    );
                    let decorator_line: String = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                        mutation.chrom, first_coord, last_coord, mut_name, 0, strand_lit,
                        first_coord, last_coord, mutation.color, 1, 1, 0,
                        tr_coords, GLYPH, mutation.color, shape
                    );
                    if let Err(e) = output_file.write(decorator_line.as_bytes()) {
                        eprintln!("Failed to write the line: {}", e);
                    }
                };
            };
        }
    }
}