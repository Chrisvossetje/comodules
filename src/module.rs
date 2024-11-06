pub trait Module {

}

pub trait ModuleMorphism {
    fn get_domain(&self) -> impl Module;
    fn get_codomain(&self) -> impl Module;

    fn get_cokernel(&self) -> impl ModuleMorphism;
}