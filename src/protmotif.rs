use ::*;

#[derive(Serialize)]
pub struct ProtMotif {
    pub seq_ct: usize,
    pub scores: Array2<f32>,
    /// sum of "worst" base at each position
    pub min_score: f32,
    /// sum of "best" base at each position
    pub max_score: f32,
}

impl ProtMotif {
    pub fn from_seqs_with_pseudocts(seqs: Vec<Vec<u8>>, pseudos: &[f32; 20]) -> ProtMotif {
        if seqs.len() == 0 {
            return ProtMotif {
                seq_ct: 0,
                scores: Array2::zeros((0, 0)),
                min_score: 0.0,
                max_score: 0.0,
            };
        }

        let seqlen = seqs[0].len();
        let mut counts = Array2::zeros((seqlen, 20));
        for i in 0..seqlen {
            for base in 0..20 {
                counts[[i, base]] = pseudos[base];
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
        let mut m = ProtMotif {
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
            for base_i in 0..20 {
                tot += self.scores[[i, base_i]];
            }
            for base_i in 0..20 {
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
            let min_sc = (0..20).map(|b| self.scores[[i, b]]).fold(
                f32::INFINITY,
                f32::min,
            );
            self.min_score += min_sc;
        }

        // score corresponding to "best" base at each position
        self.max_score = 0.0;
        for i in 0..pwm_len {
            let max_sc = (0..20).map(|b| self.scores[[i, b]]).fold(
                f32::NEG_INFINITY,
                f32::max,
            );
            self.max_score += max_sc;
        }
    }

}

impl Motif for ProtMotif {
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
        4,
        3,
        5,
        13,
        7,
        8,
        9,
        255,
        11,
        10,
        12,
        2,
        255,
        14,
        6,
        1,
        15,
        16,
        255,
        19,
        17,
        255,
        18,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        0,
        255,
        4,
        3,
        5,
        13,
        7,
        8,
        9,
        255,
        11,
        10,
        12,
        2,
        255,
        14,
        6,
        1,
        15,
        16,
        255,
        19,
        17,
        255,
        18,
        255,
        255,
        255,
        255,
        255,
    ];
    const MONOS: &'static [u8] = b"ARNDCEQGHILKMFPSTWYV";

    fn rev_lk(idx: usize) -> u8 {
        if idx >= Self::MONOS.len() {
            INVALID_MONO
        } else {
            Self::MONOS[idx]
        }
    }

    fn len(&self) -> usize {
        self.scores.dim().0
    }
    fn get_scores(&self) -> &Array2<f32> {
        &self.scores
    }
    fn get_min_score(&self) -> f32 {
        self.min_score
    }
    fn get_max_score(&self) -> f32 {
        self.max_score
    }
    fn get_bits() -> f32 {
        20f32.log2()
    }
    fn degenerate_consensus(&self) -> Vec<u8> {
        let len = self.len();
        let mut res = Vec::with_capacity(len);
        for pos in 0..len {
            let mut fracs = (0..20)
                .map(|b| (self.scores[[pos, b]], b))
                .collect::<Vec<(f32, usize)>>();
            // note: reverse sort
            fracs.sort_by(|a, b| b.partial_cmp(a).unwrap());

            res.push(if fracs[0].0 > 0.5 && fracs[0].0 > 2.0 * fracs[1].0 {
                Self::MONOS[fracs[0].1]
            } else {
                b'X'
            });
        }
        res
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_info_content() {
        let pwm = ProtMotif::from_seqs_with_pseudocts(vec![b"AAAA".to_vec()], &[0.0; 20]);
        assert_eq!(pwm.info_content(), ProtMotif::get_bits() * 4.0);
    }
}
