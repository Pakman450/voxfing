use voxbirch::voxelize;
use voxbirch::read_mol2_file;
use voxbirch::write_cluster_mol_ids;
use voxbirch::birch::VoxBirch;
use voxbirch::get_recommended_info;
use voxbirch::calc_time_breakdown;

use std::path::{Path};
use nalgebra::DMatrix;
use clap::Parser;
use std::time::{Instant};
use env_logger::{Builder};
use std::io::Write;
use log::LevelFilter;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the MOL2 file (required)
    #[arg(short, long, required = true)]
    path: String,

    /// Dimensions of the voxel grid (x, y, z), comma separated
    #[arg(short, long, value_delimiter = ',', default_value = "20,20,20")]
    dims: Vec<usize>,

    /// Resolution of the voxel grid in Angstroms
    #[arg(short, long, default_value_t = 2.0)]
    resolution: f32,

    /// Target origin x0 y0 z0 via comma separated string
    #[arg(short, long, value_delimiter = ',', default_value = "0.0,0.0,0.0")]
    origin: Vec<f32>,

    /// Threshold of similarity
    #[arg(short, long, default_value_t = 0.65)]
    threshold: f32,

    /// Number of max branches
    #[arg(short, long, default_value_t = 50)]
    max_branches: usize,

    /// Clustered mol ids output name
    #[arg(long, default_value = "./clustered_mol_ids.txt")]
    output_file_path: String,

    /// Do not condense voxel grids. Leaving this out condenses grids.
    #[arg(long)]
    no_condense: bool,

    /// Verbosity level. -v means level 1, -vvv means level 3 
    #[arg(short, long, action = clap::ArgAction::Count, default_value_t = 0)]
    verbosity: u8,

}

fn init_logging(verbosity: u8) {
    let mut builder = Builder::new();

    match verbosity {
        0 => builder.filter_level(LevelFilter::Warn),   // default
        1 => builder.filter_level(LevelFilter::Info),
        2 => builder.filter_level(LevelFilter::Debug),
        _ => builder.filter_level(LevelFilter::Trace),
    };

    builder.format(|buf, record| {

        let level_style = buf.default_level_style(record.level());
        let level = level_style.value(record.level());

        let file = record.file().unwrap_or("unknown");
        let line = record.line().unwrap_or(0);

        writeln!(
            buf,
            "[{} {}:{} {}] {}",
            level,
            file,
            line,
            record.target(),
            record.args()
        )
    });
    
    builder.init();
}

fn main() {

    // Get the current time
    let start_time = Instant::now();

    // Print the ASCII art
    voxbirch::ascii::print_ascii_art();
        
    // Argument unpacking
    let args = Args::parse();
    let file_path= args.path;
    let dimx = args.dims[0];
    let dimy = args.dims[1];
    let dimz = args.dims[2];
    let x0 = args.origin[0];
    let y0 = args.origin[1];       
    let z0 = args.origin[2];
    let resolution = args.resolution;
    let threshold = args.threshold;
    let max_branches = args.max_branches;
    let clustered_mol_id_string = args.output_file_path;
    let no_condense = args.no_condense;
    init_logging(args.verbosity);

    // Read MOL2 file
    let path = Path::new(&file_path);
    let l_mols = read_mol2_file(path).expect("Failed to read MOL2 file");

    // Print some input info
    println!("################################################");
    println!("MOL2 file path: {}", file_path);
    println!("Number of molecules read: {}", l_mols.len());
    println!("Voxel Grid Dimensions: {} x {} x {}", dimx, dimy, dimz);
    println!("Voxel Grid Resolution: {}", resolution);
    println!("Voxel Grid Origin: ({}, {}, {})", x0, y0, z0);
    println!("Thresold: {}", threshold);
    println!("Max branches: {}", max_branches);
    println!("Condense Voxel Grids: {}", !no_condense);
    println!("################################################");

    println!("\nGrabbing recommended voxelization parameters...");
    // Give user the recommended origin values for placing voxels.
    let (
        min_x, 
        min_y, 
        min_z, 
        need_x, 
        need_y, 
        need_z,
        need_x_user,
        need_y_user,
        need_z_user
    ) = get_recommended_info(&l_mols, resolution, x0, y0, z0);
    println!("The recommended origin: {},{},{}", min_x.floor(), min_y.floor(), min_z.floor());

    println!(
        "Minimal voxel grid dimensions to cover all molecules\n\t(from absolute origin {}): {},{},{}",
        format!("{:.3},{:.3},{:.3}", min_x, min_y, min_z), need_x, need_y, need_z
    );
    println!(
        "Required dims from provided origin ({:.3},{:.3},{:.3}): {},{},{}", 
        x0, y0, z0, 
        need_x_user, need_y_user, need_z_user
    );

    println!("\nVoxelizing...");
    // Voxelization of molecules's xyz's
    let grids = voxelize(&l_mols, [dimx, dimy, dimz], resolution, x0, y0, z0); 

    let voxelize_duration: std::time::Duration = start_time.elapsed();

    // Get the number of rows (which is the number of VoxelGrids)
    let num_rows = grids.len();

    // Get the number of cols (which is the number of voxels)
    let num_cols;

    if !no_condense {
        num_cols = grids[0].condensed_data.len(); 
    } else{
        num_cols = grids[0].data.len(); 
    }
    
    println!("Shape of data: ({} molecules, {} voxels)", num_rows, num_cols);

    // Create the DMatrix with the correct size
    let mut input_matrix: DMatrix<f32> = DMatrix::zeros(
        num_rows, 
        num_cols
    );

    if !no_condense {
        for (i, grids) in grids.iter().enumerate() {
            for (j, &value) in grids.condensed_data.iter().enumerate() {
                input_matrix[(i, j)] = value as f32; // Convert u8 to f32 and assign
            }
        }
    } else{
        for (i, grids) in grids.iter().enumerate() {
            for (j, &value) in grids.data.iter().enumerate() {
                input_matrix[(i, j)] = value as f32; // Convert u8 to f32 and assign
            }
        }
    }

    let mut vb = VoxBirch::new(
        threshold, 
        max_branches
    );

    if !input_matrix.iter().any(|&x| x != 0.0) {
        panic!("All of your voxels for all rows have 0.0s ");
    }

    // start clustering
    vb.fit(
        &input_matrix, 
        grids.iter().map(|g| g.title.clone()).collect()
    );

    // Get results after clustering. 
    let cluster_mol_ids: Vec<Vec<String>> = vb.get_cluster_mol_ids();
    let num_clusters = cluster_mol_ids.len();
    let path_cluster_ids: String = String::from(clustered_mol_id_string);
    let write_to_path = Path::new(&path_cluster_ids);
    let _ = write_cluster_mol_ids(&write_to_path, &cluster_mol_ids);

    // Get the breakdown of elapsed time
    let (
        vox_hours,
        vox_minutes,
        vox_seconds,
        vox_milliseconds
    ) = calc_time_breakdown(&voxelize_duration);

    let total_duration: std::time::Duration = start_time.elapsed();
    let (
        tot_hours,
        tot_minutes,
        tot_seconds,
        tot_milliseconds
    ) = calc_time_breakdown(&total_duration);

    println!(
"\nFinished
Summary statistics:
Total number of clusters: {}
Elapsed time for Voxelization: {}h {}m {}s {}ms
Elapsed time for Clustering: {}h {}m {}s {}ms
Total Elapsed time: {}h {}m {}s {}ms",
        num_clusters,
        
        vox_hours, 
        vox_minutes,
        vox_seconds, 
        vox_milliseconds,
        
        tot_hours - vox_hours, 
        tot_minutes - vox_minutes, 
        tot_seconds - vox_seconds, 
        tot_milliseconds - vox_milliseconds,
        
        tot_hours, 
        tot_minutes, 
        tot_seconds, 
        tot_milliseconds
    );

}
