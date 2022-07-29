use clap::Parser;
use faimm::IndexedFasta;
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

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
    if let Ok(lines) = read_lines(args.vcf) {
        for line in lines {
            lineno = lineno + 1;
            if let Ok(ip) = line {
                if ip.chars().next().unwrap() == '#' {
                    continue;
                }
                let s = ip.to_string();
                let (chr_name, pos, ref_allele) = get_columns(&s);
                let chr_index = fa.fai().tid(chr_name).expect("Cannot find chr in index");
                let v = fa
                    .view(chr_index, pos, pos + ref_allele.len())
                    .expect("Cannot get .fa view");

                let needle = v.to_string().to_lowercase();
                let haystack = ref_allele.to_lowercase();
                if needle != haystack {
                    println!("Failed line {}: {}", lineno, ip);
                    err = 1;
                }
            }
        }
    }

    if err != 0 {
        Err("Found some mismatching lines")
    } else {
        println!("No mismatching lines found");
        Ok(())
    }
}

// from https://doc.rust-lang.org/rust-by-example/std_misc/file/read_lines.html#read_lines
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<flate2::read::GzDecoder<File>>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;

    let reader = io::BufReader::new(GzDecoder::new(file));
    Ok(reader.lines())
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
