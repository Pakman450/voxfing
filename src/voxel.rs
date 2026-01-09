use crate::file_io::VoxMol;

#[derive(Clone)]
pub struct VoxelGrid {
    pub title : String,
    pub dims: [usize; 3],
    pub data: Vec<u8>,    
}

impl VoxelGrid {
    pub fn new(dims: [usize; 3]) -> Self {
        let size = dims[0] * dims[1] * dims[2];
        Self {
            title: String::new(),
            dims,
            data: vec![0; size],
        }
    }

    // This function computes the 1D index in the voxel grid from 3D coordinates
    pub fn voxel_index(&self, x: usize, y: usize, z: usize) -> usize {
        x + self.dims[0] * (y + self.dims[1] * z)
    }

}

// Example "voxelizer": fills grid with a single value
// FIX: currently only for positive numbers
// I have to Shift coordinates to positive first
// I should find the minimium coordinates first and shift accordingly
pub fn voxelize(
    l_mols: Vec<VoxMol>,
    dims: [usize; 3],
    rs: f32,
    x0: f32,
    y0: f32,
    z0: f32
) -> Vec::<VoxelGrid> {

    // create a list of voxel grids based on number of molecules

    let mut grids = Vec::<VoxelGrid>::new();

    for mol in l_mols.iter() {

        let mut grid = VoxelGrid::new(dims);

        grid.title = mol.title.clone();

        for atom_idx in 0..mol.num_atoms() {
            let x = mol.x[atom_idx];
            let y = mol.y[atom_idx];
            let z = mol.z[atom_idx];


            let ix = ((x - x0) / rs).floor() as usize;
            let iy = ((y - y0) / rs).floor() as usize;
            let iz = ((z - z0) / rs).floor() as usize;
        
            let index = grid.voxel_index(ix, iy, iz);

            grid.data[index] += 1;

        }

        grids.push(grid);
    }

    let grid = &grids[0];  // for now just return the first grid
    // print out the number of occupied voxels
    let occupied_voxels = grid.data.iter().filter(|&&v| v > 0).count();
    println!("Number of occupied voxels: {}", occupied_voxels);

    // print out the sum of all voxel values
    let voxel_sum: usize = grid.data.iter().map(|&v| v as usize).sum();
    println!("Sum of all voxel values: {}", voxel_sum);

    grids
}


pub fn get_recommended_info(l_mols: &Vec<VoxMol>, resolution: f32, x0: f32, y0: f32, z0: f32) -> (
    f32, f32, f32, usize, usize, usize, usize, usize, usize
) {
    let mut l_vals = Vec::<f32>::new();

    for mol in l_mols {
        l_vals.extend(mol.x.clone())
    }

    let min_x = l_vals.iter().cloned().filter(|&x| !x.is_nan()).min_by(|a, b| a.partial_cmp(b).unwrap()).expect("No valid minimum found (possibly due to NaN values).");

    l_vals.clear();

    for mol in l_mols {
        l_vals.extend(mol.y.clone())
    }

    let min_y = l_vals.iter().cloned().filter(|&x| !x.is_nan()).min_by(|a, b| a.partial_cmp(b).unwrap()).expect("No valid minimum found (possibly due to NaN values).");

    l_vals.clear();
    for mol in l_mols {
        l_vals.extend(mol.z.clone())
    }
    let min_z = l_vals.iter().cloned().filter(|&x| !x.is_nan()).min_by(|a, b| a.partial_cmp(b).unwrap()).expect("No valid minimum found (possibly due to NaN values).");
        // Compute maximum coordinates across all molecules so we can compute
    // the minimal dimensions required to cover them from the recommended origin
    let mut xs: Vec<f32> = Vec::new();
    let mut ys: Vec<f32> = Vec::new();
    let mut zs: Vec<f32> = Vec::new();

    for mol in l_mols {
        xs.extend(mol.x.iter().cloned());
        ys.extend(mol.y.iter().cloned());
        zs.extend(mol.z.iter().cloned());
    }

    let max_x = xs.iter().cloned().filter(|v| !v.is_nan()).max_by(|a, b| 
        a.partial_cmp(b).unwrap()).expect("No valid maximum found (possibly due to NaN values)."
    );
    let max_y = ys.iter().cloned().filter(|v| !v.is_nan()).max_by(|a, b| 
        a.partial_cmp(b).unwrap()).expect("No valid maximum found (possibly due to NaN values)."
    );
    let max_z = zs.iter().cloned().filter(|v| !v.is_nan()).max_by(|a, b| 
        a.partial_cmp(b).unwrap()).expect("No valid maximum found (possibly due to NaN values)."
    );

    // span and required voxels (inclusive of boundary). Use floor(span / resolution) + 1
    let span_x = (max_x - min_x).max(0.0);
    let span_y = (max_y - min_y).max(0.0);
    let span_z = (max_z - min_z).max(0.0);

    let need_x = (span_x / resolution).floor() as usize + 1;
    let need_y = (span_y / resolution).floor() as usize + 1;
    let need_z = (span_z / resolution).floor() as usize + 1;

    // Also compute required dims if using the user-provided origin (x0,y0,z0)
    let span_x_user = (max_x - x0).max(0.0);
    let span_y_user = (max_y - y0).max(0.0);
    let span_z_user = (max_z - z0).max(0.0);

    let need_x_user = (span_x_user / resolution).floor() as usize + 1;
    let need_y_user = (span_y_user / resolution).floor() as usize + 1;
    let need_z_user = (span_z_user / resolution).floor() as usize + 1;

    (
        min_x, min_y, min_z, // minimum coordinates
        need_x, need_y, need_z, // required dims from absolute origin
        need_x_user, need_y_user, need_z_user // required dims from user-provided origin
    )
}



