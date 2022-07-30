use clap::Parser;
use faimm::IndexedFasta;
use rust_htslib::bgzf::Reader;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, Read};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    #[clap(short, long, value_parser)]
    fasta: String,

    #[clap(short, long, value_parser)]
    vcf: String,
}

fn main() -> Result<(), &'static str> {
    let args = Args::parse();
    let fa = IndexedFasta::from_file(args.fasta).expect("Error opening fasta file");
    let mut lineno = 0;
    let mut err = 0;
    let lines = get_reader(&args.vcf.to_string());

    for line in lines {
        lineno = lineno + 1;
        if let Ok(ip) = line {
            if ip.chars().next().unwrap() == '#' {
                continue;
            }
            let s = ip.to_string();
            let (chr_name, pos, ref_allele) = get_columns(&s);
            let chr_index = fa.fai().tid(chr_name).expect("Cannot find chr in index");
            let fasta = fa
                .view(chr_index, pos, pos + ref_allele.len())
                .expect("Cannot get .fa view")
                .to_string();

            if fasta.to_lowercase() != ref_allele.to_lowercase() {
                println!("Failed line {}: {}", lineno, ip);
                err = 1;
            }
        }
    }
    println!("Lines processed: {}", lineno);

    if err != 0 {
        Err("Found some mismatching lines")
    } else {
        println!("No mismatching lines found");
        Ok(())
    }
}

fn get_reader(path: &str) -> io::Lines<BufReader<Box<dyn Read>>> {
    let file_type = path
        .split(".")
        .collect::<Vec<&str>>()
        .last()
        .unwrap()
        .clone();

    match file_type {
        "gz" => {
            let mm: Box<dyn Read> = Box::new(Reader::from_path(path).unwrap());
            BufReader::new(mm).lines()
        }
        _t => {
            let file = File::open(path).unwrap();
            let mm: Box<dyn Read> = Box::new(file);
            BufReader::new(mm).lines()
        }
    }
}

// get chr_name, pos, and ref_allele avoiding splitting the whole line (could be many genotypes)
fn get_columns(s: &str) -> (&str, usize, &str) {
    let mut iter = s.match_indices('\t');
    let j = iter.nth(0).unwrap().0;
    let chr_name = &s[0..j];
    let k = iter.nth(0).unwrap().0;
    let posstr = &s[j + 1..k];

    let i = iter.nth(0).unwrap().0 + 1;
    let l = iter.nth(0).unwrap().0;
    let ref_allele = &s[i..l];
    let pos: usize = posstr.parse::<usize>().unwrap() - 1;
    (chr_name, pos, ref_allele)
}
