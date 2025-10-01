use bigtools::BigBedRead;
use clap::Parser;
use colorconv::Color as _Color;
use fxhash::FxHashMap;
use serde_json;
use std::cmp::{max, min};
use std::io::{BufRead, Write};
use std::fs::File;
use std::path::Path;
use std::str::FromStr;

use rust::read;

use svgdom::{
    AttributeId, AttributeValue, Color, Document, 
    ElementId, Length, LengthList, LengthUnit, 
    Node, NodeType, ViewBox
};
// use svg::write;

// const STDOUT: &str = "stdout";
const FREE_SPACE: f64 = 1.0;

const LONG_DESCRIPTION: &str = "";

struct Coords{
    start: u32,
    end: u32
}

#[derive(Debug)]
struct Dimensions{
    height: f64,
    width: f64
}

struct Plot {
    doc: Document,
    head_node: Node,
    dimensions: Dimensions,
    height_mod: f64,
    counter: u8
}

// #[derive(Deserialize, Debug)]
// struct Plot {

// }

fn recursive_addition<'a>(
    // from: &mut Document, 
    to: &mut Document,
    donor_node: &'a mut Node,
    recipient_node: &'a mut Node,
    dims: &mut Dimensions,
    mut height_mod: f64,
    plot_counter: u8,
    species: &str,
    add_species: bool
) {
    let mut is_svg_node: bool = false;
    for mut child in donor_node.children() {
        if child.node_type() == NodeType::Element {
            is_svg_node |= child.tag_id().unwrap() == ElementId::Svg;
            let mut new_node = to.create_element(
                child.tag_id().unwrap().clone()
            );
            // iterate over the elements
            for attr in  child.attributes().iter() {
                // ViewBox attribute; in TOGA2, projection-wise plots are guaranteed 
                // to have only one ViewBox 
                if let AttributeValue::ViewBox(view) = attr.value {
                    let box_height = view.h;
                    if plot_counter > 0 {dims.height += FREE_SPACE}
                    height_mod = dims.height;
                    dims.height += box_height;
                    let box_width = view.w;
                    dims.width = dims.width.max(box_width);
                    continue
                }
                if child.is_root() {
                    continue
                }
                if plot_counter > 1 {
                    let attr_id = attr.id().unwrap();
                    if attr_id == AttributeId::Y || attr_id == AttributeId::Y1 || attr_id == AttributeId::Y2 {
                        let adj_len = match &attr.value {
                            AttributeValue::Length(x) => {
                                AttributeValue::Length(Length::new(x.num + height_mod, LengthUnit::None))
                            },
                            AttributeValue::LengthList(x) => {
                                AttributeValue::LengthList(
                                    LengthList(
                                        x
                                        .iter()
                                        .map(|f| Length::new(f.num + height_mod, LengthUnit::None))
                                        .collect::<Vec::<Length>>()
                                    )
                                )
                            }
                            _ => {
                                continue
                            }
                        };
                        new_node.set_attribute((attr_id, adj_len));
                        continue
                    }
                }
                new_node.set_attribute(attr.clone());
            }
            if child.has_children() {
                recursive_addition(
                    to, 
                    &mut child,
                    if is_svg_node {recipient_node}  else {&mut new_node}, 
                    dims, 
                    height_mod, 
                    plot_counter,
                    species,
                    add_species
                );
            }
            if !is_svg_node {recipient_node.append(new_node.clone())};
        } else if child.node_type() == NodeType::Text {
            let mut upd_text = child.text().clone();
            if add_species && upd_text.contains('#') {
                upd_text = format!("{} {}", species, upd_text);
            }
            let new_node = to.create_node(NodeType::Text,upd_text);
            recipient_node.append(new_node.clone());
        }
    }
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
/// Extracts SVG plots for given projections from a TOGA2 BigBed file
struct Args{
    /// Input Plot specification in JSON format
    #[arg()]
    json: String,

    /// Path to output directory. The plots will be saved to this directory 
    /// according to the input JSON file
    #[arg()]
    output: String,

    /// If set, species names are not added to the projection IDs in the plot
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    do_not_add_species_names: bool,

    /// Background fill color; if not set, outputs the combined plot with transparent background
    #[arg(long, short = 'c')]
    fill_color: Option<String>,
}

fn main() {
    let args = Args::parse();

    // check if the output directory exists
    match Path::new(&args.output).try_exists() {
        Ok(true) => {},
        Ok(false) => {
            panic!("Output directory \"{}\" does not exist", args.output)
        },
        Err(_) => {
            panic!("Cannot confirm the existence of the output directory \"{}\"", args.output)
        }
    }

    // read input json
    let input_data: serde_json::Value = serde_json::from_reader(read(args.json))
        .expect("Failed to read the JSON file");
    // bigBed2projections mapper
    let mut path2projs: FxHashMap<String, Vec<String>> = FxHashMap::default();
    // bed2projections mapper
    let mut bed2projs: FxHashMap<String, Vec<String>> = FxHashMap::default();
    // proj2plot mapper
    let mut proj2plots: FxHashMap<String, Vec<String>> = FxHashMap::default();
    // proj2path mapper
    let mut proj2path: FxHashMap<String, String> = FxHashMap::default();
    // species names
    let mut proj2species: FxHashMap<String, String> = FxHashMap::default();
    // output plots
    let mut plots: FxHashMap<String, Plot> = FxHashMap::default();
    // first level iteration: individual plots
    for (plot_name, plot) in input_data.as_object().unwrap() {
        // second level: individual output directories
        for (species, species_data) in plot.as_object().unwrap() {
            // bigBed file path
            let path = species_data
                .get("bigbed")
                .expect(
                    &format!("Missing \"bigbed\" field for species \"{}\", plot \"{}\"", species, plot_name)
                )
                .as_str()
                .expect(
                    &format!("Improper \"bigbed\" field formatting for species \"{}\", plot \"{}\"", species, plot_name)
                )
                .to_string();
            // bed file path
            let bed_file = species_data
                .get("bed")
                .expect(
                    &format!("Missing \"bed\" field  for species \"{}\", plot \"{}\""
                    , species, plot)
                )
                .as_str()
                .expect(
                    &format!(
                        "Improper \"bed\" field formatting for species \"{}\", plot \"{}\"", species, plot
                    )
                )
                .to_string();
            let projections = species_data
                .get("projections")
                .expect(
                    &format!("Missing \"projections\" field for species \"{}\", plot \"{}\""
                    , species, plot_name)
                )
                .as_array()
                .expect(
                    &format!("\"projections\" field for species \"{}\" (plot \"{}\") is not a valid array"
                    , species, plot_name)
                );
            // parse the Bed files as well
            for projection in projections {
                path2projs
                    .entry(path.clone())
                    .and_modify(|x| 
                        x.push(projection.as_str().unwrap().to_string()))
                    .or_insert(vec![projection.as_str().unwrap().to_string()]);
                bed2projs
                    .entry(bed_file.clone())
                    .and_modify(
                        |x| x.push(projection.as_str().unwrap().to_string())
                    )
                    .or_insert(vec![projection.as_str().unwrap().to_string()]);
                proj2plots
                    .entry(projection.as_str().unwrap().to_string())
                    .and_modify(|x| x.push(plot_name.clone()))
                    .or_insert(vec![plot_name.clone()]);
                if !proj2path.contains_key(projection.as_str().unwrap()) {
                    proj2path.insert(projection.as_str().unwrap().to_string(), path.clone());
                }
                if proj2species.contains_key(projection.as_str().unwrap()) {
                    panic!("Projection {} is attributed to more than one species!", projection.as_str().unwrap());
                }
                proj2species.insert(projection.as_str().unwrap().to_string(), species.clone());
            }   
        }
        plots.insert(
            plot_name.clone(),
            {
                let mut combined_plot = Document::new();
                let svg = combined_plot.create_element(ElementId::Svg);
                combined_plot.root().append(svg.clone());
                Plot {
                    doc: combined_plot,
                    head_node: svg,
                    dimensions: Dimensions {height: 0.0, width: 0.0},
                    height_mod: 0.0,
                    counter: 1
                }
            }
        );
    }

    // step 2: retrieve the coordinates from the Bed files
    let mut bigbed2chrom2coords: FxHashMap<String, FxHashMap<String, Coords>> = FxHashMap::default();
    // TODO: chrom2coords should be tied to a specific BigBed file
    // meanwhile, we also need a 
    for (bed, projections) in bed2projs {
        for (i, line_) in read(bed).lines().enumerate() {
            if let Ok(line) = line_ {
                // let bed_record = parse_bed(line, 12, true);
                let proj_comps: Vec<&str> = line.split('\t').collect::<Vec<&str>>();
                if proj_comps.len() < 4 {
                    panic!(
                        "Improper Bed file formatting at line {}; expected at least 4 columns, got {}",
                        i+1, proj_comps.len()
                    )
                }
                let proj_name = proj_comps[3];
                // if !projections.contains(proj_name) {continue}
                if !projections.iter().any(|x| x == proj_name) {continue}
                let chrom: String = String::from(proj_comps[0]);
                let start: u32 = proj_comps[1]
                    .parse::<u32>()
                    .expect(
                        &format!(
                            "Improper Bed file formatting at line {}: {} is not a valid 64-bit integer", 
                            i+1, proj_comps[1]
                        )
                    );
                let end: u32 = proj_comps[2]
                    .parse::<u32>()
                    .expect(
                        &format!(
                            "Improper Bed file formatting at line {}: {} is not a valid 64-bit integer", 
                            i+1, proj_comps[2]
                        )
                    );
                let bigbed = proj2path
                    .get(proj_name)
                    .expect(&format!("Undefined BigBed file for projection {}", proj_name));
                bigbed2chrom2coords
                    .entry(bigbed.clone())
                    .and_modify(|x|{
                        x.entry(chrom.clone())
                            .and_modify(|y| {
                                y.start = min(y.start, start);
                                y.end = max(y.end, end);
                            })
                            .or_insert(Coords { start: start, end: end });
                    })
                    .or_insert(
                        FxHashMap::from_iter([(chrom.clone(), Coords {start, end})])
                    );
            }
        }
    }

    // step 3: Extract the plots from the BigBed files
    for (bigbed, proj_data) in bigbed2chrom2coords {
        let projections = path2projs
            .get(&bigbed)
            .expect(&format!("BigBed file {} does not contain any assigned projections", bigbed));
        if let Ok(mut bb) = BigBedRead::open_file(bigbed) {
            // iterate over chromosome intervals
            for (chrom, coords) in proj_data {
                if let Ok(entry_inter) = bb.get_interval(&chrom, coords.start, coords.end) {
                    // find the entries corresponding to the sought projections
                    for entry_ in entry_inter {
                        if let Ok(entry) = entry_ {
                            let rest_comps = entry.rest.split('\t').collect::<Vec<&str>>();
                            let name: &str = rest_comps[0];
                            if !projections.contains(&name.to_string()) {continue}
                            // the projection was requested for plotting
                            let plot = rest_comps[27];
                            // parse the SVG string as a document
                            let tree = Document::from_str(plot).unwrap();
                            // get the species name
                            let species = proj2species
                                .get(name)
                                .expect(&format!("No species attribution found for projection {}", name));
                            // add the projection data to each of the plots it belongs to
                            for plot_name in proj2plots.get(name).unwrap() {
                                let plot_data = plots
                                    .get_mut(plot_name)
                                    .expect(
                                        &format!("No plot stub found for plot \"{}\"", plot_name)
                                    );
                                recursive_addition(
                                    &mut plot_data.doc, 
                                    &mut tree.root(),//svg_element().unwrap(),//tree.root(), 
                                    &mut plot_data.head_node, 
                                    &mut plot_data.dimensions, 
                                    plot_data.height_mod, 
                                    plot_data.counter,
                                    species,
                                    !args.do_not_add_species_names
                                );
                                plot_data.counter += 1;
                            }
                        }
                    }
                }
            }
        }
    }

    for (plot_name, plot) in plots {
        plot.doc.svg_element().unwrap().set_attribute(
            (
                AttributeId::ViewBox, 
                AttributeValue::ViewBox(ViewBox::new(0.0, 0.0, plot.dimensions.width, plot.dimensions.height))
            )
        );
        plot.doc.svg_element().unwrap().set_attribute(
            (
                AttributeId::Height,
                AttributeValue::Length(Length::new(plot.dimensions.height, LengthUnit::None)),
            )
        );
        plot.doc.svg_element().unwrap().set_attribute(
            (
                AttributeId::Width,
                AttributeValue::Length(Length::new(plot.dimensions.width, LengthUnit::None)),
            )
        );
        match args.fill_color {
            Some(ref color) => {
                let rgb = _Color::from_str(&color)
                    .expect(&format!("Invalid background color provided: {}", &color))
                    .rgb;
                plot.doc.svg_element().unwrap().set_attribute(
                    (
                        AttributeId::Fill,
                        AttributeValue::Color(Color::new(rgb[0], rgb[1], rgb[2]))
                    )
                );
            },
            None => {}
        };
        println!("plot={}", plot.doc.to_string());
        let output_path = Path::new("").join(&args.output).join(format!("{plot_name}.svg"));
        let mut output_handle = Box::new(File::create(&output_path).unwrap()) as Box<dyn Write>;
        if let Err(e) = output_handle.write(plot.doc.to_string().as_bytes()) {
            eprintln!("Failed to write the line: {}", e);
        }
    }

    // save the results
    // let mut output_file = match args.output.as_str() {
    //     "stdout" => {Box::new(stdout()) as Box<dyn Write>},
    //     _ => {
    //         let path = Path::new(&args.output);
    //         Box::new(File::create(&path).unwrap()) as Box<dyn Write>
    //     }
    // };
    // if let Err(e) = output_file.write(combined_plot.to_string().as_bytes()) {
    //     eprintln!("Failed to write the line: {}", e);
    // }
}