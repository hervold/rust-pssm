use ::*;

/// monomers, ie, DNA bases
static MONOMERS: &'static [u8] = b"ATGC";
/// default pseudocount - used to prevent 0 tallies
const DEF_PSEUDO: f32 = 0.5;


struct DNAMotif {
    pub seq_ct: usize,
    pub scores: Array2<f32>,
    /// sum of "worst" base at each position
    pub min_score: f32,
    /// sum of "best" base at each position
    pub max_score: f32,
}

impl DNAMotif {
    pub fn from_seqs_with_pseudocts(seqs: Vec<Vec<u8>>, pseudos: &[f32; 4]) -> DNAMotif {
        // null case
        if seqs.len() == 0 {
            return DNAMotif {
                seq_ct: 0,
                scores: Array2::zeros((0, 0)),
                min_score: 0.0,
                max_score: 0.0,
            };
        }

        let seqlen = seqs[0].len();
        let mut counts = Array2::zeros((seqlen, 4));
        for i in 0..seqlen {
            for base in 0..4 {
                counts[[i,base]] = pseudos[base];
            }
        }

        for seq in seqs.iter() {
            if seq.len() != seqlen {
                panic!("inconsistent sequence lengths");
            }

            for (idx, base) in seq.iter().enumerate() {
                counts[[idx, Self::lookup(*base).expect("DNA base")]] += 1.0;
            }
        }
        let mut m = DNAMotif {
            seq_ct: seqs.len(),
            scores: counts,
            min_score: 0.0,
            max_score: 0.0,
        };
        m.normalize();
        m.calc_minmax();
        m
    }

    // helper function -- normalize self.scores
    fn normalize(&mut self) {
        for i in 0..self.len() {
            let mut tot: f32 = 0.0;
            // FIXME: slices would be cleaner
            for base_i in 0..4 {
                tot += self.scores[[i, base_i]];
            }
            for base_i in 0..4 {
                self.scores[[i, base_i]] /= tot;
            }
        }
    }

    // helper function
    fn calc_minmax(&mut self) {
        let pwm_len = self.len();

        // score corresponding to sum of "worst" bases at each position
        // FIXME: iter ...
        self.min_score = 0.0;
        for i in 0..pwm_len {
            // can't use the regular min/max on f32, so we use f32::min
            let min_sc = (0..4).map(|b| self.scores[[i, b]]).fold(
                f32::INFINITY,
                f32::min,
            );
            self.min_score += min_sc;
        }

        // score corresponding to "best" base at each position
        self.max_score = 0.0;
        for i in 0..pwm_len {
            let max_sc = (0..4).map(|b| self.scores[[i, b]]).fold(
                f32::NEG_INFINITY,
                f32::max,
            );
            self.max_score += max_sc;
        }
    }
}

impl Motif for DNAMotif {
    const LK: [u8; 127] = [
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        0,
        255,
        3,
        255,
        255,
        255,
        2,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        1,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        0,
        255,
        3,
        255,
        255,
        255,
        2,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        1,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
    ];

    fn len(&self) -> usize {
        self.scores.dim().0
    }

    fn get_monos(&self) -> &[u8] {
        MONOMERS
    }

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
        for start in 0..seq.len() - pwm_len + 1 {
            let m: Vec<f32> = (0..pwm_len)
                .map(|i| {
                    self.scores[[i, Self::lookup(seq[start + i]).expect("raw lookup")]]
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
        if self.max_score == self.min_score {
            return None;
        }

        let (best_start, best_score, best_m) = self.raw_score(&seq);

        if best_score < 0.0 || best_score - self.min_score < 0.0 ||
            self.max_score - self.min_score < 0.0
        {
            println!(
                "@@ found problem: best = {}; num = {} - {} = {}; demon = {}",
                best_score,
                best_score,
                self.min_score,
                best_score - self.min_score,
                self.max_score - self.min_score
            );
        }
        Some(ScoredPos {
            loc: best_start,
            sum: (best_score - self.min_score) / (self.max_score - self.min_score),
            scores: best_m,
        })
    }

    /// derived from
    /// https://github.com/biopython/biopython/blob/master/Bio/motifs/matrix.py#L205
    fn degenerate_consensus(&self) -> Vec<u8> {
        fn two(_a: u8, _b: u8) -> u8 {
            let (a, b) = if _b > _a { (_a, _b) } else { (_b, _a) };
            match (a, b) {
                (b'A', b'C') => b'M',
                (b'A', b'G') => b'R',
                (b'A', b'T') => b'W',
                (b'C', b'G') => b'S',
                (b'C', b'T') => b'Y',
                (b'G', b'T') => b'K',
                _ => panic!("unrecognized bases: {}, {}", a, b),
            }
        }
        let len = self.len();
        let mut res = Vec::with_capacity(len);
        for pos in 0..len {
            let mut fracs = (0..4)
                .map(|b| (self.scores[[pos, b]], b))
                .collect::<Vec<(f32, usize)>>();
            // note: reverse sort
            fracs.sort_by(|a, b| b.partial_cmp(a).unwrap());

            res.push(if fracs[0].0 > 0.5 && fracs[0].0 > 2.0 * fracs[1].0 {
                MONOMERS[fracs[0].1]
            } else if 4.0 * (fracs[0].0 + fracs[1].0) > 3.0 {
                two(MONOMERS[fracs[0].1], MONOMERS[fracs[1].1])
            } else if fracs[3].0 < EPSILON {
                match MONOMERS[fracs[3].1] {
                    b'T' => b'V',
                    b'G' => b'H',
                    b'C' => b'D',
                    b'A' => b'B',
                    _ => panic!("weirdo triplet"),
                }
            } else {
                b'N'
            });
        }
        res
    }
}

/// calculate scores matrix from a list of equal-length sequences
/// use DEF_PSEUDO as default pseudocount
impl From<Vec<Vec<u8>>> for DNAMotif {
    fn from(seqs: Vec<Vec<u8>>) -> Self {
        DNAMotif::from_seqs_with_pseudocts(seqs, &[DEF_PSEUDO, DEF_PSEUDO, DEF_PSEUDO, DEF_PSEUDO])
    }
}

impl From<Array2<f32>> for DNAMotif {
    fn from(scores: Array2<f32>) -> Self {
        let mut m = DNAMotif {
            seq_ct: 0,
            scores: scores,
            min_score: 0.0,
            max_score: 0.0,
        };
        m.normalize();
        m.calc_minmax();
        m
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn simple_pwm() {
        let pwm = DNAMotif::from(vec![
            b"AAAA".to_vec(),
            b"TTTT".to_vec(),
            b"GGGG".to_vec(),
            b"CCCC".to_vec(),
        ]);
        assert_eq!(pwm.scores, Array2::from_elem((4, 4), 0.25));
    }
    #[test]
    fn find_motif() {
        let pwm = DNAMotif::from(vec![b"ATGC".to_vec()]);
        let seq = b"GGGGATGCGGGG";
        if let Some(ScoredPos { ref loc, ref sum, .. }) = pwm.score(seq) {
            assert_eq!(*loc, 4);
            assert_eq!(*sum, 1.0);
        } else {
            assert!(false);
        }
    }
}
