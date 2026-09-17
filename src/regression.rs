//!Bit-exact state dumps used to prove a refactor did not move any number.
//!Inert unless main.rs sets a dump directory.

use crate::pipes::InteriorMethod;
use std::collections::BTreeMap;
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::path::PathBuf;

///Creates the dump directory and clears the trace left by any previous run.
pub fn init(dir: &str) {
    std::fs::create_dir_all(dir).expect("could not create dump directory");
    File::create(trace_path(dir)).expect("could not create trace file");
}

///Iterations that get a full state snapshot, chosen to localize a first divergence.
pub fn is_checkpoint(it: u32) -> bool {
    matches!(it, 0 | 1 | 10 | 100)
}

///Writes every real cell of every pipe as raw f64 bits, one line per cell.
/// Bits rather than decimals so the comparison is exact.
pub fn dump_state(dir: &str, it: u32, pipes: &BTreeMap<usize, InteriorMethod>) {
    let mut out = String::new();

    for pipe in pipes.values() {
        let s = pipe.solver().state();
        for c in 0..s.n_real {
            let col = s.q1.column(s.first + c);
            out.push_str(&format!(
                "p{} c{} {:016x} {:016x} {:016x}\n",
                pipe.id(),
                c,
                col[0].to_bits(),
                col[1].to_bits(),
                col[2].to_bits()
            ));
        }
    }

    let path = PathBuf::from(dir).join(format!("it{it:06}.state"));
    File::create(&path)
        .and_then(|mut f| f.write_all(out.as_bytes()))
        .unwrap_or_else(|e| panic!("could not write {}: {e}", path.display()));
}

///Appends the raw bits of the current time and timestep for one iteration.
/// A matching trace proves the loops line up before any state is compared.
pub fn trace_step(dir: &str, it: u32, t: f64, dt: f64) {
    let line = format!("{it} {:016x} {:016x}\n", t.to_bits(), dt.to_bits());
    let path = trace_path(dir);
    OpenOptions::new()
        .append(true)
        .open(&path)
        .and_then(|mut f| f.write_all(line.as_bytes()))
        .unwrap_or_else(|e| panic!("could not append {}: {e}", path.display()));
}

fn trace_path(dir: &str) -> PathBuf {
    PathBuf::from(dir).join("steps.trace")
}

