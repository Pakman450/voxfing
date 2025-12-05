use crate::file_io::VoxMol;


pub struct VoxelGrid {
    pub dims: [usize; 3],
    pub data: Vec<u8>,    
}

impl VoxelGrid {
    pub fn new(dims: [usize; 3]) -> Self {
        let size = dims[0] * dims[1] * dims[2];
        Self {
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
pub fn voxelize(
    l_mols: Vec<VoxMol>,
    dims: [usize; 3],
    rs: f32,
    x0: f32,
    y0: f32,
    z0: f32,
) -> VoxelGrid {

    let mut grid = VoxelGrid::new(dims); // 4×4×4 grid


    for mol in l_mols.iter() {
        for atom_idx in 0..mol.num_atoms() {
            let x = mol.x[atom_idx];
            let y = mol.y[atom_idx];
            let z = mol.z[atom_idx];

            let ix = ((x - x0) / rs).floor() as usize;
            let iy = ((y - y0) / rs).floor() as usize;
            let iz = ((z - z0) / rs).floor() as usize;
            let index = grid.voxel_index(ix, iy, iz);
            grid.data[index] = 1; // Mark voxel as occupied
        }
    }

    // print out the number of occupied voxels
    let occupied_voxels = grid.data.iter().filter(|&&v| v > 0).count();
    println!("Number of occupied voxels: {}", occupied_voxels);

    grid
}




