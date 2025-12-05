use voxelizer::voxelize;
use voxelizer::read_mol2_file;

use std::env;
use std::path::{Path};


fn main() {

    let args: Vec<String> = env::args().collect(); // get CLI args as Vec


    let file_path: String = args[1].parse().unwrap_or_else(|_| {
        eprintln!("Invalid value: {}", args[1]);
        std::process::exit(1);
    });
    if args.len() < 2 {
        eprintln!("Usage: voxelize <value>");
        return;
    }
    
    let path = Path::new(&file_path);

    match read_mol2_file(path) {
        Ok(l_mols) => {
            for mol in l_mols.iter() {
                println!("{:?}", mol.x);
                println!("{:?}", mol.y);
                println!("{:?}", mol.z);
                println!("Number of atoms: {:?}", mol.num_atoms());
            }
        }
        Err(e) => eprintln!("Error reading file: {}", e),
    }

    let grid = voxelize(7); // fill with value 7
    println!("Grid dims: {:?}", grid.dims);
    println!("First 8 values: {:?}", &grid.data[..8]);
}
