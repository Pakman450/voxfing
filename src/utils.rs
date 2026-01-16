
pub fn calc_time_breakdown (duration_mark: &std::time::Duration) -> (
    u64,
    u64,
    u64,
    u32
) {

    let total_secs = duration_mark.as_secs();
    let hours = total_secs / 3600;
    let minutes = (total_secs % 3600) / 60;
    let seconds = total_secs % 60;
    let milliseconds = duration_mark.subsec_millis();

    (
        hours,
        minutes,
        seconds,
        milliseconds
    )
}