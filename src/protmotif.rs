use ::*;

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
    const MONOS: &'static [u8] = b"ARNDCEQGHILKMFPSTWYV";

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
        vec![]
    }
}
