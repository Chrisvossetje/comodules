pub type AGeneratorIndex = u16;
pub type ComoduleIndexType = u32;
pub type CoalgebraIndexType = u16;

pub type ComoduleIndex<G> = (G, ComoduleIndexType);
pub type CoalgebraIndex<G> = (G, CoalgebraIndexType);

pub type CofreeBasis<G> = (CoalgebraIndex<G>, AGeneratorIndex);

pub type UniGradingType = u8;
