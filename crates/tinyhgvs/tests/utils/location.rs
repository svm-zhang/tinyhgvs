use tinyhgvs::{Interval, Location, NucleotideCoordinate, ProteinCoordinate};

pub trait LocationExt<T> {
    fn known_start(&self) -> &T;
    fn left_bp(&self) -> &Interval<T>;
    fn right_bp(&self) -> &Interval<T>;
}

impl LocationExt<NucleotideCoordinate> for Location<NucleotideCoordinate> {
    fn known_start(&self) -> &NucleotideCoordinate {
        self.start()
            .expect("expect start coordinate for known location")
    }

    fn left_bp(&self) -> &Interval<NucleotideCoordinate> {
        self.l_interval()
            .expect("expect uncertain left breakpoint interval")
    }

    fn right_bp(&self) -> &Interval<NucleotideCoordinate> {
        self.r_interval()
            .expect("expect uncertain right breakpoint interval")
    }
}

impl LocationExt<ProteinCoordinate> for Location<ProteinCoordinate> {
    fn known_start(&self) -> &ProteinCoordinate {
        self.start()
            .expect("expect start coordinate for known location")
    }

    fn left_bp(&self) -> &Interval<ProteinCoordinate> {
        self.l_interval()
            .expect("expect uncertain left breakpoint interval")
    }

    fn right_bp(&self) -> &Interval<ProteinCoordinate> {
        self.r_interval()
            .expect("expect uncertain right breakpoint interval")
    }
}

pub trait IntervalExt<T> {
    fn interval_end(&self) -> &T;
}

impl<T> IntervalExt<T> for Interval<T> {
    fn interval_end(&self) -> &T {
        self.end.as_ref().expect("expect interval end")
    }
}

pub fn assert_five_prime_intron_coord(coord: &NucleotideCoordinate) {
    assert!(coord.is_intronic());
    assert!(coord.is_cds_start_anchored());
    assert!(!coord.is_cds_end_anchored());
    assert!(!coord.is_five_prime_utr());
    assert!(!coord.is_three_prime_utr());
}

pub fn assert_three_prime_intron_coord(coord: &NucleotideCoordinate) {
    assert!(coord.is_intronic());
    assert!(!coord.is_cds_start_anchored());
    assert!(coord.is_cds_end_anchored());
    assert!(!coord.is_five_prime_utr());
    assert!(!coord.is_three_prime_utr());
}

pub fn assert_five_prime_utr_coord(coord: &NucleotideCoordinate) {
    assert!(!coord.is_intronic());
    assert!(coord.is_cds_start_anchored());
    assert!(!coord.is_cds_end_anchored());
    assert!(coord.is_five_prime_utr());
    assert!(!coord.is_three_prime_utr());
}

pub fn assert_three_prime_utr_coord(coord: &NucleotideCoordinate) {
    assert!(!coord.is_intronic());
    assert!(!coord.is_cds_start_anchored());
    assert!(coord.is_cds_end_anchored());
    assert!(!coord.is_five_prime_utr());
    assert!(coord.is_three_prime_utr());
}
