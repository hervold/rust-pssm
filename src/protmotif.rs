use ::*;

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
    const width: usize = 20;
    const lk: [u8; 127] = [255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
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

}
