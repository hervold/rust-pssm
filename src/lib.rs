
#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate serde_derive;
extern crate bio;
extern crate ndarray;

use bio::utils::IntoTextIterator;
use std::f32;
use std::convert::From;
use ndarray::prelude::Array2;

mod dnamotif;
mod protmotif;

pub use dnamotif::DNAMotif;
pub use protmotif::ProtMotif;

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

pub trait Motif {
    /// lookup table mapping monomer -> index
    const LK: [u8; 127] = [INVALID_MONO; 127];
    const MONOS: &'static [u8] = b"";

    // FIXME: this is a hack to work around some trouble I had with generics
    fn from_scores(scores: Array2<f32>) -> Self;

    /// use lk to find index; enforce boundaries
    fn lookup(mono: u8) -> Option<usize> {
        if mono >= 127 {
            None
        } else {
            let idx = Self::LK[mono as usize];
            if idx == INVALID_MONO {
                None
            } else {
                Some(idx as usize)
            }
        }
    }
    /// reverse lookup: given index, return monomer (eg, 0 -> A)
    fn rev_lk(idx: usize) -> u8;

    fn len(&self) -> usize;

    /// represent highly conserved bases
    fn degenerate_consensus(&self) -> Vec<u8>;
    fn get_scores(&self) -> &Array2<f32>;
    /// sum of "worst" base at each position
    fn get_min_score(&self) -> f32;
    /// sum of "best" base at each position
    fn get_max_score(&self) -> f32;
    fn get_bits() -> f32;


    // standard PSSM scoring is calibrated to (1) pattern len and (2) the min and
    //     max possible scores
    // this is just a dumb sum of matching bases
    fn raw_score<'a, T: IntoTextIterator<'a>>(&self, seq_it: T) -> (usize, f32, Vec<f32>) {
        let pwm_len = self.len();

        let mut best_start = 0;
        let mut best_score = -1.0;
        let mut best_m = Vec::new();
        // we have to look at slices, so a simple iterator won't do
        let seq = seq_it.into_iter().cloned().collect::<Vec<u8>>();
        let scores = self.get_scores();
        for start in 0..seq.len() - pwm_len + 1 {
            let m: Vec<f32> = (0..pwm_len)
                .map(|i| {
                    scores[[i, Self::lookup(seq[start + i]).expect("raw lookup")]]
                })
                .collect();
            let tot = m.iter().sum();
            if tot > best_score {
                best_score = tot;
                best_start = start;
                best_m = m;
            }
        }
        (best_start, best_score, best_m)
    }

    /// apply PSM to sequence, finding the offset with the highest score
    /// return None if sequence is too short
    /// see:
    ///   MATCHTM: a tool for searching transcription factor binding sites in DNA sequences
    ///   Nucleic Acids Res. 2003 Jul 1; 31(13): 3576–3579
    ///   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC169193/
    ///
    fn score<'a, T: IntoTextIterator<'a>>(&self, seq_it: T) -> Option<ScoredPos> {
        let pwm_len = self.len();
        let seq = seq_it.into_iter().cloned().collect::<Vec<u8>>();
        if seq.len() < pwm_len {
            return None;
        }
        let min_score = self.get_min_score();
        let max_score = self.get_max_score();

        if max_score == min_score {
            return None;
        }

        let (best_start, best_score, best_m) = self.raw_score(&seq);

        Some(ScoredPos {
            loc: best_start,
            sum: (best_score - min_score) / (max_score - min_score),
            scores: best_m,
        })
    }

    /// roughly the inverse of Shannon Entropy
    /// adapted from the information content described here:
    ///    https://en.wikipedia.org/wiki/Sequence_logo#Logo_creation
    fn info_content(&self) -> f32 {
        fn ent<'a, I>(probs: I) -> f32
        where
            I: Iterator<Item = &'a f32>,
        {
            probs
                .map(|p| if *p == 0.0 {
                    0.0
                } else {
                    -1.0 * *p * p.log(2.0)
                })
                .sum()
        }
        let bits = Self::get_bits();
        let scores = self.get_scores();
        let mut tot = 0.0;
        for row in scores.genrows() {
            tot += bits - ent(row.iter());
        }
        tot
    }
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
