
#[macro_use]
extern crate serde_derive;
extern crate bio;
extern crate ndarray;

use bio::utils::{IntoTextIterator, TextIterator};
use std::f32;
use std::convert::From;
use std::default;
use ndarray::prelude::{Array, Array2, Axis};

mod dnamotif;
//mod protmotif;

pub const EPSILON: f32 = 1e-5;
pub const INVALID_MONO: u8 = 255;

/// represents motif score
#[derive(Debug, Clone, Serialize)]
pub struct ScoredPos {
    pub loc: usize,
    pub sum: f32,
    pub scores: Vec<f32>,
}

impl Default for ScoredPos {
    fn default() -> ScoredPos {
        ScoredPos {
            loc: 0,
            sum: f32::NEG_INFINITY,
            scores: Vec::new(),
        }
    }
}

trait Motif {
    /// lookup table mapping monomer -> index
    const lk: [u8; 127] = [INVALID_MONO; 127];
    /// use lk to find index; enforce boundaries
    fn lookup(mono: u8) -> Option<usize> {
        if mono >= 127 {
            None
        } else {
            let idx = Self::lk[mono as usize];
            if idx == INVALID_MONO { None } else { Some(idx as usize) }
        }
    }

    fn len(&self) -> usize;
    fn get_monos(&self) -> &[u8];
    /// represent highly conserved bases
    fn degenerate_consensus(&self) -> Vec<u8>;

    /// simply the sum of matching monomers; can serve as a dumber info_content
    fn raw_score<'a, T: IntoTextIterator<'a>>(&self, seq: T) -> (usize, f32, Vec<f32>);

    /// apply PSM to sequence, finding the offset with the highest score
    /// return None if sequence is too short
    /// see:
    ///   MATCHTM: a tool for searching transcription factor binding sites in DNA sequences
    ///   Nucleic Acids Res. 2003 Jul 1; 31(13): 3576–3579
    ///   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC169193/
    ///
    fn score<'a, T: IntoTextIterator<'a>>(&self, seq: T) -> Option<ScoredPos>;
    /*
    /// roughly the inverse of Shannon Entropy
    /// adapted from the information content described here:
    ///    https://en.wikipedia.org/wiki/Sequence_logo#Logo_creation
    fn info_content(&self) -> f32;*/
}


    /*
    /// update scores field based on counts & pseudocounts
    pub fn calc_scores(&mut self) {
        match self.seqs {
            Some(_) => {
                let seqs = self.seqs.as_ref().expect("calc_scores called without seqs");
                for seq_i in 0..seqs[0].len() {
                    let mut tot: f32 = 0.0;
                    // FIXME: slices should be cleaner
                    for base_i in 0..4 {
                        tot += self.counts.as_ref().expect(
                            "calc_scores called without counts",
                        )
                            [[seq_i, base_i]] as f32;
                        tot += self.pseudoct[base_i];
                    }
                    let mut scores = self.scores.subview_mut(Axis(0), seq_i);
                    for base_i in 0..4 {
                        scores[base_i] =
                            (self.counts.as_ref().expect(
                                "calc_scores called without counts",
                            )
                                 [[seq_i, base_i]] as f32 +
                                 self.pseudoct[base_i]) / tot;
                    }
                }
            }
            _ => (),
        }
        self.calc_minmax();
    }
     */
