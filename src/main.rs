use clap::Parser;
use faimm::IndexedFasta;
use rust_htslib::bgzf::Reader;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, Read};
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
    let (lineno, err) = verify(&args.fasta, &args.vcf);
    println!("Lines processed: {}", lineno);

    if err != 0 {
        Err("Found some mismatching lines")
    } else {
        println!("No mismatching lines found");
        Ok(())
    }
}

fn verify(fasta: &str, vcf: &str) -> (usize, usize) {
    let fa = IndexedFasta::from_file(fasta).expect("Error opening fasta file");
    let mut lineno = 0;
    let mut err = 0;
    let lines = get_reader(vcf);

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
                println!("Failed line {}: {} {}", lineno, ip, fasta);
                err = 1;
            }
        }
    }
    return (lineno, err);
}

fn get_reader(path: &str) -> io::Lines<BufReader<Box<dyn Read>>> {
    match Path::new(&path).extension() {
        None => {
            let file = File::open(path).unwrap();
            let mm: Box<dyn Read> = Box::new(file);
            BufReader::new(mm).lines()
        }
        Some(os_str) => match os_str.to_str() {
            Some("gz") => {
                let mm: Box<dyn Read> = Box::new(Reader::from_path(path).unwrap());
                BufReader::new(mm).lines()
            }
            _t => {
                let file = File::open(path).unwrap();
                let mm: Box<dyn Read> = Box::new(file);
                BufReader::new(mm).lines()
            }
        },
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plaintext_vcf() {
        let (lineno, err) = verify("test/volvox.fa", "test/volvox.filtered.vcf");
        assert_eq!(lineno, 73);
        assert_eq!(err, 0);
    }

    #[test]
    fn gzip_vcf() {
        let (lineno, err) = verify("test/volvox.fa", "test/volvox.filtered.vcf.gz");
        assert_eq!(lineno, 73);
        assert_eq!(err, 0);
    }

    #[test]
    fn bgzip_vcf() {
        let (lineno, err) = verify("test/volvox.fa", "test/volvox.filtered.vcf.gz");
        assert_eq!(lineno, 73);
        assert_eq!(err, 0);
    }

    #[test]
    #[should_panic]
    fn missing_vcf() {
        verify("test/volvox.fa", "missing.gz");
    }

    #[test]
    #[should_panic]
    fn missing_fasta() {
        verify("missing.fa", "missing.gz");
    }
}
