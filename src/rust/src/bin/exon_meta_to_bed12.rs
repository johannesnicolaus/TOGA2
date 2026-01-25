use clap::Parser;
use cubiculum::extract::extract::to_line;
use cubiculum::structs::structs::{BedEntry, Named};
use fxhash::FxHashMap;
use phf::phf_map;
use std::fs::File;
use std::io::{BufRead, stdout, Write};
use std::path::Path;
use rust::read;

#[derive(Parser, Debug)]
#[command(version = "1.0", about, long_about = None)]
/// Reconstructs the BED12 file from the exon metadata file, ignoring non-Present exons.
/// WARNING: Currently ignores the fragmented projections
struct Args {
    /// Input exon metadata file; can be gzip-compressed. Not providing a value or setting value 
    /// to "stdin" implies reading data from standard input 
    #[arg(long, short = 'i', default_value_t = String::from("stdin"))]
    input: String,

    /// Loss summary file; necessary for loss status extraction. 
    /// Not providing a value or setting value 
    /// to "stdin" implies reading data from standard input 
    #[arg(long, short = 'l', default_value_t = String::from("stdin"))]
    loss_file: String,

    /// Output file; if a value is not provided or set to "stdout", the 
    /// result will be written to standard output
    #[arg(long, short = 'o', default_value_t = String::from("stdout"))]
    output: String
}

const EXON_META_HEADER: &str = "projection";
const INTACT: &str = "I";
const PROJECTION: &str = "PROJECTION";

const LOSS2COLOR: phf::Map<&'static str, &'static str> = phf_map!{
    "PP" => "250,50,200",
    "PG" => "159,129,112",
    "L" => "255,50,50",
    "M" => "130,130,130",
    "UL" => "255,160,120",
    "PI" => "0,200,255",
    "I" => "0,0,200",
    "FI" => "0,0,100",
};

fn main() {
    let args = Args::parse();
    let mut proj2loss: FxHashMap<String, String> = FxHashMap::default();
    for (i, line_) in read(args.loss_file).lines().enumerate() {
        if let Ok(line) = line_ {
            let fields: Vec<&str> = line.split('\t').collect::<Vec<&str>>();
            if fields.len() == 0 {continue}
            if fields.len() != 3 {
                panic!(
                    "Improper loss file formatting at line {}; expected 3 fields, got {}",
                    i + 1, fields.len()
                )
            }
            if fields[0] != PROJECTION {continue}
            proj2loss.insert(String::from(fields[1]), String::from(fields[2]));
        }
    }
    let mut bed_records: FxHashMap<String, BedEntry> = FxHashMap::default();
    let mut fragmented_bed_records: FxHashMap<String, FxHashMap<String, BedEntry>> = FxHashMap::default();
    let mut exon2chain: FxHashMap<String, FxHashMap<String, Vec<u16>>> = FxHashMap::default();
    for (i, line_) in read(args.input).lines().enumerate() {
        if let Ok(line) = line_ {
            let fields: Vec<&str> = line.split('\t').collect::<Vec<&str>>();
            if fields.len() == 0 {continue}
            if fields.len() < 22 {
                panic!(
                    "Improper formatting at exon meta file line {}; expected at least 22 fields, got {}", 
                    i + 1, fields.len()
                )
            }
            if fields[0] == EXON_META_HEADER {continue}
            let projection: &str = fields[0];
            let exon: u16 = fields[1]
                .parse::<u16>()
                .expect(
                    &format!(
                        "Improper formatting at exon meta line {}: Exon number value is not a valid numeric value",
                        i+1
                    )
                );
            let chain: String = String::from(fields[2]);
            let exon_status: &str = fields[7];
            if exon_status != INTACT {continue}
            let chrom: &str = fields[3];
            let start: u64 = fields[4]
                .parse::<u64>()
                .expect(
                    &format!(
                        "Improper formatting at exon meta line {}: Exon start value is not a valid numeric value",
                        i + 1
                    )
                );
            let end: u64 = fields[5]
            .parse::<u64>()
            .expect(
                &format!(
                    "Improper formatting at exon meta line {}: Exon end value is not a valid numeric value",
                    i + 1
                )
            );
            if end - start == 0 {continue}
            let loss_status = proj2loss
                .get(projection)
                .expect(
                    &format!(
                        "Projection {} is missing from the loss file", projection
                    )
                );
            let rgb: String = String::from(*LOSS2COLOR.get(loss_status).unwrap());
            let strand: bool = fields[6] == "+";
            let exon_entry = BedEntry::bed12(
                String::from(chrom), 
                start, 
                end, 
                String::from(projection), 
                String::from("1000"), 
                strand,
                start,
                end,
                rgb,
                1,
                vec![end - start],
                vec![0]
            );
            // insert/graft the resulting exon
            if projection.contains(',') {
                // assign the exon to the respective chain fragment
                // if it's the first time the projection is encountered, create the storage HashMap
                if let None = exon2chain.get(projection) {
                    fragmented_bed_records.insert(
                        String::from(projection), 
                        FxHashMap::from_iter(vec![(chain.clone(), exon_entry)])
                    );
                    exon2chain.insert(
                        String::from(projection), 
                        FxHashMap::from_iter(vec![(chain.clone(), vec![exon])])
                    );
                } else {
                    // insert and/or graft the value
                    fragmented_bed_records
                        .entry(String::from(projection))
                        .and_modify(|y|{
                            y
                            .entry(chain.clone())
                            .and_modify(|x|{
                                x.graft(
                                    exon_entry.clone(),
                                    true,
                                    true,
                                    true,
                                    true,
                                    false,
                                    false
                                );
                            })
                            .or_insert(exon_entry);
                        });
                    exon2chain
                        .entry(String::from(projection))
                        .and_modify(|y| {
                            y
                                .entry(chain)
                                .and_modify(|x| x.push(exon))
                                .or_insert(vec![exon]);
                        });
                }
            } else {
                bed_records
                    .entry(String::from(projection))
                    .and_modify(|x| {
                        x.graft(
                            exon_entry.clone(), 
                            true, 
                            true, 
                            true, 
                            true, 
                            false, 
                            false
                        );
                    })
                    .or_insert(exon_entry);
            }
        }
    }

    let mut output_file = match args.output.as_str() {
        "stdout" => {Box::new(stdout()) as Box<dyn Write>},
        _ => {
            let path = Path::new(&args.output);
            Box::new(File::create(&path).unwrap()) as Box<dyn Write>
        }
    };
    for (projection, bed_entry) in bed_records {
        let result_line = to_line(&bed_entry, 12).expect(
            &format!("Failed to convert entry for {} into a BED line", projection)
        );
        if let Err(e) =  writeln!(output_file, "{}", result_line) {
            eprintln!("Failed to write the line: {}", e);
        }
    }
    for (projection, chains) in fragmented_bed_records {
        let mut chains_sorted: Vec<String> = chains
            .keys()
            .cloned()
            .collect();
        chains_sorted.sort_by(|a, b| {
            exon2chain.get(&projection).unwrap().get(a).unwrap().iter().min().unwrap().cmp(
                exon2chain.get(&projection).unwrap().get(b).unwrap().iter().min().unwrap()
            ) 
        });
        let mut chain_counter: u8 = 1;
        for chain in chains_sorted {
            let mut entry = chains.get(&chain).unwrap().clone();
            entry.update_name(&format!("{}${}", entry.name().unwrap(), chain_counter));
            let result_line = to_line(&entry, 12).expect(
                &format!("Failed to convert entry for {} into a BED line", projection)
            );
            if let Err(e) =  writeln!(output_file, "{}", result_line) {
                eprintln!("Failed to write the line: {}", e);
            }
            chain_counter += 1;
        }
    }
}