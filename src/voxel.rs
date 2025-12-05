

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
}



// Example "voxelizer": fills grid with a single value
pub fn voxelize(value: u8) -> VoxelGrid {

    let mut grid = VoxelGrid::new([4, 4, 4]); // 4×4×4 grid
    grid.data.fill(value);
    grid
}




