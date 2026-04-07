/// FASTA file reading with gzip support

use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use flate2::read::GzDecoder;
use anyhow::Result;

/// Represents a sequence record from FASTA
#[derive(Debug, Clone)]
pub struct FastaRecord {
    pub name: String,
    pub sequence: String,
}

/// Reader abstraction for both plain and gzipped FASTA
pub enum FastaReader {
    Plain(BufReader<File>),
    Gzipped(BufReader<GzDecoder<File>>),
}

impl FastaReader {
    /// Create a new FASTA reader, detecting gzip format from file extension
    pub fn new(path: &str) -> Result<Self> {
        let file = File::open(path)?;
        
        if path.ends_with(".gz") {
            let gz = GzDecoder::new(file);
            Ok(FastaReader::Gzipped(BufReader::new(gz)))
        } else {
            Ok(FastaReader::Plain(BufReader::new(file)))
        }
    }

    /// Iterate over all records in the FASTA file
    pub fn records(&mut self) -> FastaIterator {
        FastaIterator {
            reader: self,
            buffer: String::new(),
            current_name: String::new(),
            current_seq: String::new(),
            finished: false,
        }
    }
}

/// Iterator for FASTA records
pub struct FastaIterator<'a> {
    reader: &'a mut FastaReader,
    buffer: String,
    current_name: String,
    current_seq: String,
    finished: bool,
}

impl<'a> Iterator for FastaIterator<'a> {
    type Item = Result<FastaRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        loop {
            self.buffer.clear();
            
            // Read line
            let bytes_read = match self.reader {
                FastaReader::Plain(ref mut r) => {
                    match r.read_line(&mut self.buffer) {
                        Ok(n) => n,
                        Err(e) => return Some(Err(anyhow::Error::from(e))),
                    }
                }
                FastaReader::Gzipped(ref mut r) => {
                    match r.read_line(&mut self.buffer) {
                        Ok(n) => n,
                        Err(e) => return Some(Err(anyhow::Error::from(e))),
                    }
                }
            };

            if bytes_read == 0 {
                // End of file
                self.finished = true;
                if !self.current_name.is_empty() {
                    return Some(Ok(FastaRecord {
                        name: self.current_name.clone(),
                        sequence: self.current_seq.clone(),
                    }));
                }
                return None;
            }

            let line = self.buffer.trim_end();

            if line.starts_with('>') {
                // Header line
                if !self.current_name.is_empty() {
                    // Return previous record
                    let name = self.current_name.clone();
                    let seq = self.current_seq.clone();
                    
                    self.current_name = line[1..].to_string();
                    self.current_seq.clear();
                    
                    return Some(Ok(FastaRecord {
                        name,
                        sequence: seq,
                    }));
                } else {
                    // First header
                    self.current_name = line[1..].to_string();
                }
            } else if !line.is_empty() {
                // Sequence line
                self.current_seq.push_str(line);
            }
        }
    }
}

/// Load all FASTA sequences into memory
pub fn load_fasta(path: &str) -> Result<Vec<FastaRecord>> {
    let mut records = Vec::new();
    let mut reader = FastaReader::new(path)?;

    for record_result in reader.records() {
        let record = record_result?;
        records.push(record);
    }

    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fasta_parsing() {
        // TODO: Add real FASTA parsing tests
    }
}
