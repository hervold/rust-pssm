use ::*;


lazy_static! {
    static ref LOG20: f32 = 20f32.log2();
}
static AMINO_ACIDS: &'static [u8] = b"ARNDCEQGHILKMFPSTWYV";


struct ProtMotif {
    pub seqs: Option<Vec<Vec<u8>>>,
    pub counts: Option<Array2<u16>>,
    pub pseudoct: [f32; 20],
    pub seq_ct: usize,
    pub scores: Array2<f32>,
    /// sum of "worst" base at each position
    pub min_score: f32,
    /// sum of "best" base at each position
    pub max_score: f32,
}

impl Motif for ProtMotif {
    const LK: [u8; 127] = [255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                           255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                           255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                           255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                           255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                           255, 255, 255, 255, 255, 0, 255, 4, 3, 5, 13, 7,
                           8, 9, 255, 11, 10, 12, 2, 255, 14, 6, 1, 15,
                           16, 255, 19, 17, 255, 18, 255, 255, 255, 255, 255, 255,
                           255, 0, 255, 4, 3, 5, 13, 7, 8, 9, 255, 11,
                           10, 12, 2, 255, 14, 6, 1, 15, 16, 255, 19, 17,
                           255, 18, 255, 255, 255, 255, 255];
    fn get_monos(&self) -> &[u8] {
        AMINO_ACIDS
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
    fn degenerate_consensus(&self) -> Vec<u8> {
        vec![]
    }
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

        let mut tot = 0.0;
        for row in self.scores.genrows() {
            tot += *LOG20 - ent(row.iter());
        }
        tot
    }

}
