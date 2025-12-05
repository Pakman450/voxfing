use voxelizer::voxelize;
use voxelizer::read_mol2_file;

use std::env;
use std::path::{Path};


fn main() {

    let args: Vec<String> = env::args().collect(); // get CLI args as Vec

    // if args.len() < 2 {
    //     eprintln!("Usage: voxelize <value>");
    //     return;
    // }

    let file_path: String = args[1].parse().unwrap_or_else(|_| {
        eprintln!("Invalid value: {}", args[1]);
        std::process::exit(1);
    });

    let dimx = args[2].parse().unwrap_or_else(|_| {
        eprintln!("Invalid value: {}", args[2]);
        std::process::exit(1);
    });
    let dimy = args[3].parse().unwrap_or_else(|_| {
        eprintln!("Invalid value: {}", args[3]);
        std::process::exit(1);
    });
    let dimz = args[4].parse().unwrap_or_else(|_| {
        eprintln!("Invalid value: {}", args[4]);
        std::process::exit(1);
    });


    let resolution: f32 = args[5].parse().unwrap_or_else(|_| {
        eprintln!("Invalid value: {}", args[5]);
        std::process::exit(1);
    });


    
    let path = Path::new(&file_path);

    // match read_mol2_file(path) {
    //     Ok(l_mols) => {
    //         for mol in l_mols.iter() {
    //             println!("{:?}", mol.x);
    //             println!("{:?}", mol.y);
    //             println!("{:?}", mol.z);
    //             println!("Number of atoms: {:?}", mol.num_atoms());
    //         }
    //     }
    //     Err(e) => eprintln!("Error reading file: {}", e),
    // }

    let l_mols = read_mol2_file(path).expect("Failed to read MOL2 file");



    let grid = voxelize(l_mols, [dimx, dimy, dimz], resolution, 0.0, 0.0, 0.0); 
    println!("First 8 values: {:?}", &grid.data[..8]);
}
