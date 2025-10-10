use anyhow::{Context, Error, Result};
use chaintools::cmap::chain::Chain;
use cubiculum::merge::merge::intersection;
use cubiculum::structs::structs::{Coordinates, Interval, Named};
use flate2::read::MultiGzDecoder;
use std::cmp::{max, min, Ord};
use std::fmt::Debug;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Main reading function; reads either a plain text file or a gzipped file line-by-line
pub fn read<T>(file: T) -> Box<dyn BufRead> 
    where 
    T: AsRef<Path>
{
if file.as_ref().extension().unwrap() == "gz" {
    let path = File::open(file).unwrap();
    return Box::new(BufReader::new(MultiGzDecoder::new(path))) as Box<dyn BufRead>
} else {
    let path = File::open(file).unwrap();
    return Box::new(BufReader::new(path)) as Box<dyn BufRead>
}
}

/// An auxiliary function for obtaining coding sequence boundaries in both reference and query
/// 
/// # Arguments
/// `intervals` - A Vector of Interval-like objects. All objects are expected to have defined coordinates
/// and names
/// 
/// # Returns
/// A Result object containing a Tuple of Option<Interval>, where the first interval stands for 
/// coordinates in the reference (min(first_exon_start, first_block_start), max(last_exon_end, last_block_end)) 
/// and the next one contains terminal CDS base projections in the query
/// 
pub fn get_cds_boundaries_<'a, T>(
    chain: &Chain, intervals: &'a mut Vec<T>
) -> Result<(Option<Interval>, Option<Interval>)> 
where 
    T: Coordinates + Named + Debug
{
    let mut ref_interval: Interval = Interval::new();
    let mut query_interval: Interval = Interval::new();

    // extract strand in the reference
    let ref_strand = chain.refs.strand == '+';

    // sort the input intervals; if reference strand is negative, revert the 
    match ref_strand {
        true => {
            intervals.sort_by(
                |a, b| if a.start().unwrap() == b.start().unwrap() {
                    a.end().unwrap().cmp(&b.end().unwrap())
                } else {
                    a.start().unwrap().cmp(&b.start().unwrap())
                }
            );
        },
        false => {
            intervals.sort_by(
                |a, b| if a.start().unwrap() == b.start().unwrap() {
                    a.end().unwrap().cmp(&b.end().unwrap()).reverse()
                } else {
                    a.start().unwrap().cmp(&b.start().unwrap()).reverse()
                }
            );
        }
    };
    // define the total span for the input intervals
    let mut min_start: u64 = *intervals[if ref_strand {0} else {intervals.len() - 1}]
        .start()
        .with_context(||
            {"Cannot assess coverage for intervals with undefined coordinates"}
        )?;
    // note, however,  that the elements are sorted by the start coordinate alone,
    // so the last element must not necessarily end farthest
    let mut max_end: u64 = *intervals[if ref_strand {intervals.len() - 1} else {0}]
        .end().with_context(||
            {"Cannot assess coverage for intervals with undefined coordinates"}
        )?;
    // create a smart iteration index; iteration will always start from this interval
    let mut curr: usize = 0;
    // record the current interval's end coordinate; this will ensure that the iterator will never
    // skip the nested intervals
    let mut curr_start: u64 = *intervals[0].start().with_context(||
        {"Cannot assess coverage for intervals with undefined coordinates"}
    )?;

    // create a smart iteration index; iteration will always start from this interval
    let mut curr: usize = 0;
    // record the current interval's end coordinate; this will ensure that the iterator will never
    // skip the nested intervals
    let mut curr_end: u64 = *intervals[0].end().with_context(||
        {"Cannot assess coverage for intervals with undefined coordinates"}
    )?;

    // initialize the variables standing for block coordinates
    // in this case, only the ref coordinates matter
    let mut r_start: u64;
    let r_end: u64;
    match ref_strand {
        true => {
            r_start = chain.refs.start;
            r_end = chain.refs.end;
        },
        false => {
            // NOTE: Below, r_start is defined as GREATER coordinate;
            // this is intended, since in the negative ref strand case 
            // the iterator moves upstream along the ref genome
            r_start = chain.refs.size - chain.refs.start;
            r_end = chain.refs.size - chain.refs.end;
        }
    };
    let mut r_block_start: u64 = r_start;
    let mut r_block_end: u64 = 0;

    // define query coordinates
    let q_strand: bool = chain.query.strand == '+';
    let mut q_start: u64 = match q_strand {
        true => chain.query.start,
        false => chain.query.size - chain.query.start
    };
    let query_size = chain.query.size;
    let mut q_block_start: u64 = 0;
    let mut q_block_end: u64 = 0;

    // store the indices of the marginal mapped intervals
    let mut start_in_ref: Option<u64> = None;
    let mut start_in_query: Option<u64> = None;
    let mut end_in_ref: Option<u64> = None;
    let mut end_in_query: Option<u64> = None;
    let mut first: Option<usize> = None;
    let mut last: Option<usize> = None;

    'outer: for b in chain.alignment.iter() {
        (r_block_start, r_block_end) = match ref_strand {
            true => {(r_start, r_start + b.size as u64)},
            false => {(r_start - b.size as u64, r_start)}
        };
        // println!(
        //     "chain_id={}, r_start={}, r_block_start={}, r_block_end={}, ref_strand={}, min_start={}, max_end={}, ref_start={}, ref_end={}",
        //     chain.id, r_start, r_block_start, r_block_end, ref_strand, 
        //     min_start, max_end,
        //     if ref_strand {chain.refs.start} else {chain.refs.size - chain.refs.start}, 
        //     if ref_strand {chain.refs.end} else {chain.refs.size - chain.refs.end}
        // );
        // continue if the first interval has not yet been reached
        if (ref_strand && r_block_end < min_start) || (!ref_strand && r_block_start >= max_end) {
            // don't forget to update the next block's start point
            r_start = if ref_strand {r_start + (b.size + b.dt) as u64} else {r_start - (b.size + b.dt) as u64};
            q_start = if q_strand {q_start + (b.size + b.dq) as u64} else {q_start - (b.size + b.dq) as u64};
            continue
        };
        // break the block loop if the last interval has been passed
        if (ref_strand && r_start > max_end) || (!ref_strand && r_block_end < min_start) {
            break
        };

        if q_strand {
            q_block_start = q_start;
            q_block_end = q_block_start + (b.size as u64);
        } else {
            q_block_start = q_start - (b.size as u64);
            q_block_end = q_start;
        }

        // if v {
        //     println!("Block: r_start={}, r_block_end={}, q_block_start={}, q_block_end={}, min_start={}", r_start, r_block_end, q_block_start, q_block_end, min_start);
        // }

        for (mut i, inter) in intervals[curr..].iter().enumerate() {
            i += curr;
            let inter_start: u64 = *inter.start().with_context(||
                {format!("Interval {} has an undefined start coordinate which cannot be mapped", i)}
            )?;
            let inter_end: u64 = *inter.end().with_context(||
                {format!("Interval {} has an undefined end coordinate which cannot be mapped", i)}
            )?;
            let name: &str = inter.name().with_context(||
                {"Interval is not named"}
            )?;

            // println!("inter_start={}, inter_end={}", inter_start, inter_end);

            // if v {
            //     println!("Interval: i={}, inter_start={}, inter_end={}", i, inter_start, inter_end);
            // }

            // chain block is upstream to the current interval;
            // since other are guaranteed to start at least in the same position,
            // the current loop can be safely exited
            if (ref_strand && r_block_end <= inter_start) || (!ref_strand && r_block_start >= inter_end) {
                // the pointer can be updated here, but only if the next block is guaranteed to lie further 
                // downstream to the previous interval;
                // since the chain block are sorted and do not overlap, the easiest way to prove it
                // is to check whether the current block's end does not end within the current interval group
                if (ref_strand && r_block_end > curr_end) || (!ref_strand && r_block_start < curr_start) {
                    curr = i;
                    if curr >= intervals.len() {break 'outer}
                }
                // potentially this is the farthest the intervals have reach so far
                // in terms of the end coordinate; unless this boundary is exceeded, 
                // the iteration start point will not be updated
                if inter_end >= curr_end {
                    // curr = i;
                    curr_end = inter_end;
                }
                if inter_start < curr_start {
                    curr_start = inter_start;
                }
                break
            }

            // chain block is downstream to the current interval;
            // nothing to do here, proceed to the next interval;
            if (ref_strand && r_start >= inter_end) || (!ref_strand && r_block_end <= inter_start) {
                // if inter_end == curr_end {
                //     curr += 1;
                // }
                continue
            };

            // TODO: Update for negative reference strand
            if let Some(x) = intersection(inter_start, inter_end, r_block_start, r_block_end) {
                // ignore zero-base intersection
                if x == 0 {continue}

                // the exact behavior depends on the iterator direction (=ref strand);
                // the code below could be potentially squeezed into a single if-block,
                // but for now let's keep the logic as explicit as possible
                if ref_strand {
                    // if reference strand is positive, the iterator progresses in 5'-3' direction in ref;
                    // the first exon-intersecting block yields the chain span start, span end is updated 
                    // at every other intersecting block 
                    if let None = start_in_ref {
                        // r_block_start is always a smaller coord in the tuple,
                        // therefore the ref strand is irrelevant for marginal ref coord estumation
                        start_in_ref = Some(max(r_block_start, inter_start));
                        // 
                        match r_block_start > inter_start {
                            // exon protrudes outside of the first intersected block; set it to the block marginal coordinate
                            true => {
                                start_in_query = if q_strand {Some(q_block_start)} else {Some(q_block_end)};
                            },
                            // exon starts within the block; project it through the block
                            false => {
                                let offset = inter_start - r_start;
                                start_in_query = match q_strand {
                                    true => {Some(min(q_block_start + offset, query_size))},
                                    false => {Some(q_block_end.checked_sub(offset).unwrap_or(0))}
                                }
                            }
                        };
                    }
                    end_in_ref = Some(min(r_block_end, inter_end));
                    match inter_end < r_block_end {
                        // exon starts within the block; project it through the block
                        true => {
                            let offset = r_block_end - inter_end;
                            end_in_query = match q_strand {
                                true => {Some(q_block_end.checked_sub(offset).unwrap_or(0))},
                                false => {Some(min(q_block_start + offset, query_size))}
                            }
                        },
                        // exon protrudes outside of the first intersected block; set it to the block marginal coordinate
                        false => {
                            let offset = inter_end - r_block_end;
                            // end_in_query = match q_strand {
                            //     true => {Some(min(q_block_end + offset, query_size))},
                            //     false => {Some(q_block_start.checked_sub(offset).unwrap_or(0))}
                            // };
                            end_in_query = if q_strand {Some(q_block_end)} else {Some(q_block_start)};
                        }
                    }
                } else {
                    // the iterator goes backwards (3'-5'); the first intersecting block marks the end 
                    // of the chain span in the reference, the start coordinate is further updated
                    if let None = end_in_ref {
                        end_in_ref = Some(min(r_block_end, inter_end));
                        // report the query span end
                        match r_block_end < inter_end {
                            // exon protrudes outside of the block; set the end to the marginal block coordinate
                            // invert the logic presented above for plus-strand
                            true => {
                                end_in_query = if q_strand {Some(q_block_start)} else {Some(q_block_end)};
                            },
                            false => {
                                // exon ends within the block; project it through the block
                                let offset = r_block_end - inter_end;
                                end_in_query = match q_strand {
                                    true => {Some(min(q_block_start + offset, query_size))},
                                    false => {Some(q_block_end.checked_sub(offset).unwrap_or(0))}
                                };
                            }
                        }
                    }
                    // update span start on every iteration
                    start_in_ref = Some(max(r_block_start, inter_start));
                    match r_block_start > inter_start {
                        // exon protrudes outside of the block; set the start to the marginal block coordinate
                        true => {
                            start_in_query = if q_strand {Some(q_block_end)} else {Some(q_block_start)}
                        },
                        false => {
                            // exon starts within the block => project the start
                            let offset = inter_start - r_block_start;
                            start_in_query = match q_strand {
                                true => {Some(q_block_end.checked_sub(offset).unwrap_or(0))},
                                false => {Some(min(q_block_start + offset, query_size))}
                            };
                        }
                    }
                }
            }
            curr_start = min(curr_start, inter_start);
            min_start = min(curr_start, min_start);
            curr_end = max(curr_end, inter_end);
            max_end = max(curr_end, max_end);
        }

        // if the last interval has been passed after the inner for-loop, break the outer one
        if curr >= intervals.len() {break}

        // otherwise, update the next block's start point
        r_start = if ref_strand {r_start + (b.size + b.dt) as u64} else {r_start - (b.size + b.dt) as u64};
        q_start = if q_strand {q_start + (b.size + b.dq) as u64} else {q_start - (b.size + b.dq) as u64};
        // and the first interval's start point
        match ref_strand {
            true => {min_start = *intervals[curr].start().unwrap();},
            false => {max_end = *intervals[curr].end().unwrap();}
        }

    }

    // return empty values if now exon was intersected by the coding blocks
    if let None = start_in_ref {
        return Ok((None, None))
    }

    // create an Interval object for the reference genome
    ref_interval.update_start(start_in_ref.unwrap());
    ref_interval.update_end(end_in_ref.unwrap());

    // same for query genome; 
    query_interval.update_start(min(start_in_query.unwrap(), end_in_query.unwrap()));
    query_interval.update_end(max(start_in_query.unwrap(), end_in_query.unwrap()));

    // println!("ref_interval={:#?}\nquery_interval={:#?}", ref_interval, query_interval);

    Ok((Some(ref_interval), Some(query_interval)))
}