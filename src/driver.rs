use crate::pipes::{InteriorMethod, PipeState, TimeIntegrator, apply_bc, rk_stage};
use nalgebra::Matrix3xX;
use std::collections::BTreeMap;

///Scratch buffers one pipe needs to be walked through the integrator's stages.
struct Registers {
    qn: Matrix3xX<f64>,
    stages: Vec<Matrix3xX<f64>>,
}

///Owns the time integrator and advances every pipe through its stages.
/// Registers live here rather than on the pipes so a residual can borrow one freely.
pub struct Driver {
    integrator: TimeIntegrator,
    regs: BTreeMap<usize, Registers>,
}

impl Driver {
    ///Allocates every stage register once, sized from the pipes it will advance.
    pub fn new(integrator: TimeIntegrator, pipes: &BTreeMap<usize, InteriorMethod>) -> Self {
        let regs = pipes
            .iter()
            .map(|(&id, pipe)| {
                let n_total = pipe.solver().state().n_total;
                let stages = (0..integrator.n_registers())
                    .map(|_| Matrix3xX::zeros(n_total))
                    .collect();

                (
                    id,
                    Registers {
                        qn: Matrix3xX::zeros(n_total),
                        stages,
                    },
                )
            })
            .collect();

        Self { integrator, regs }
    }

    ///Advances every pipe one full step, stage by stage.
    /// Leaves the decoded primitives stale, since q1 has moved.
    pub fn step(&mut self, dt: f64, pipes: &mut BTreeMap<usize, InteriorMethod>) {
        let stages = self.integrator.stages();

        for (id, pipe) in pipes.iter_mut() {
            let regs = self.regs.get_mut(id).expect("pipe has no stage registers");

            //copy the scalars out first so the buffers below can be split-borrowed
            let s = pipe.solver().state();
            let (first, n_real, n_ghost) = (s.first, s.n_real, s.n_ghost);
            let (left_bc, right_bc) = (s.left_bc, s.right_bc);
            let c = dt / s.dx;
            regs.qn.copy_from(&s.q1);

            for (k, &(a, b)) in stages.iter().enumerate() {
                //evaluate the residual at the state this stage starts from
                if k == 0 {
                    pipe.solver_mut().residual(&regs.qn);
                } else {
                    pipe.solver_mut().residual(&regs.stages[k - 1]);
                }

                let Registers { qn, stages: bufs } = &mut *regs;

                if k + 1 == stages.len() {
                    //last stage writes the answer straight into the pipe
                    let PipeState { q1, df, .. } = pipe.solver_mut().state_mut();
                    let src = if k == 0 { &*qn } else { &bufs[k - 1] };
                    rk_stage(q1, qn, src, df, a, b, c, first, n_real);
                    apply_bc(q1, first, n_real, n_ghost, left_bc, right_bc);
                } else {
                    let (done, rest) = bufs.split_at_mut(k);
                    let dst = &mut rest[0];
                    let src = if k == 0 { &*qn } else { &done[k - 1] };
                    rk_stage(dst, qn, src, &pipe.solver().state().df, a, b, c, first, n_real);
                    apply_bc(dst, first, n_real, n_ghost, left_bc, right_bc);
                }
            }
        }
    }
}
