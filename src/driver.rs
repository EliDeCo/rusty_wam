use crate::boundaries::{BoundaryCondition, boundary_state, euler_flux};
use crate::junctions::{Junction, ghost_state, interface_flux};
use crate::pipes::{BoundaryPair, InteriorMethod, PipeState, TimeIntegrator, apply_bc, rk_stage};
use nalgebra::{Matrix3xX, Vector3, Vector5};
use std::collections::BTreeMap;

///Scratch buffers one pipe needs to be walked through the integrator's stages.
struct Registers {
    qn: Matrix3xX<f64>,
    stages: Vec<Matrix3xX<f64>>,
}

///The same buffers for a junction, which is a single cell of five values.
struct JunctionRegisters {
    qn: Vector5<f64>,
    stages: Vec<Vector5<f64>>,
}

///One integrator stage for a junction cell: dst = a*q_n + b*(src + dt*R).
/// The junction is a control volume, so its scaling is dt/volume rather than dt/dx.
fn blend(
    q_n: &Vector5<f64>,
    src: &Vector5<f64>,
    df: &Vector5<f64>,
    a: f64,
    b: f64,
    dt_over_vol: f64,
) -> Vector5<f64> {
    debug_assert!((a + b - 1.0).abs() < 1e-15, "stage weights must sum to 1");

    a * q_n + b * (src - dt_over_vol * df)
}

///What drives one pipe end, resolved once so no stage has to search for it.
enum EndKind {
    Junction {
        junction: usize,
        ///index of this pipe's entry in the junction's own inlet list
        inlet: usize,
        normal: Vector3<f64>,
    },
    Prescribed(BoundaryCondition),
}

///One pipe end whose ghosts the driver fills every stage.
struct DrivenEnd {
    pipe: usize,
    ///whether this is the pipe's left end
    left: bool,
    ///padded index of the real cell this face touches
    cell: usize,
    ///state this end started in, which a non-reflecting boundary holds itself against
    reference: Vector3<f64>,
    area: f64,
    gamma: f64,
    kind: EndKind,
}

///Owns the time integrator and advances every pipe and junction through its stages.
/// Registers live here rather than on the objects so a residual can borrow one freely.
pub struct Driver {
    integrator: TimeIntegrator,
    pipe_regs: BTreeMap<usize, Registers>,
    junction_regs: BTreeMap<usize, JunctionRegisters>,
    ///what each junction hands its branches, rebuilt every stage
    bc: BTreeMap<usize, BoundaryPair>,
    ends: Vec<DrivenEnd>,
}

impl Driver {
    ///Allocates every stage register once, sized from the objects it will advance.
    pub fn new(
        integrator: TimeIntegrator,
        pipes: &BTreeMap<usize, InteriorMethod>,
        junctions: &BTreeMap<usize, Junction>,
    ) -> Self {
        let n_regs = integrator.n_registers();

        let pipe_regs = pipes
            .iter()
            .map(|(&id, pipe)| {
                let n_total = pipe.solver().state().n_total;
                let stages = (0..n_regs).map(|_| Matrix3xX::zeros(n_total)).collect();

                (
                    id,
                    Registers {
                        qn: Matrix3xX::zeros(n_total),
                        stages,
                    },
                )
            })
            .collect();

        let junction_regs = junctions
            .keys()
            .map(|&id| {
                (
                    id,
                    JunctionRegisters {
                        qn: Vector5::zeros(),
                        stages: vec![Vector5::zeros(); n_regs],
                    },
                )
            })
            .collect();

        let bc = pipes
            .keys()
            .map(|&id| (id, BoundaryPair::default()))
            .collect();
        let ends = resolve_ends(pipes, junctions);

        Self {
            integrator,
            pipe_regs,
            junction_regs,
            bc,
            ends,
        }
    }

    ///Rebuilds what each junction hands its branches for the stage about to run.
    /// Both sides of an interface take the same flux, which is what makes it conservative.
    fn refresh_boundary_data(&mut self, k: usize, junctions: &mut BTreeMap<usize, Junction>) {
        for pair in self.bc.values_mut() {
            *pair = BoundaryPair::default();
        }

        let Self {
            bc,
            pipe_regs,
            junction_regs,
            ends,
            ..
        } = self;

        for end in ends.iter() {
            let p_regs = &pipe_regs[&end.pipe];
            let p_state = if k == 0 {
                &p_regs.qn
            } else {
                &p_regs.stages[k - 1]
            };
            let interior = Vector3::new(
                p_state[(0, end.cell)],
                p_state[(1, end.cell)],
                p_state[(2, end.cell)],
            );

            let (ghost, flux) = match &end.kind {
                EndKind::Junction {
                    junction,
                    inlet,
                    normal,
                } => {
                    let j_regs = &junction_regs[junction];
                    let j_state = if k == 0 {
                        j_regs.qn
                    } else {
                        j_regs.stages[k - 1]
                    };
                    let (f_3d, f_1d) =
                        interface_flux(&j_state, &interior, normal, end.left, end.gamma);

                    junctions
                        .get_mut(junction)
                        .expect("end names a junction that does not exist")
                        .pipes[*inlet]
                        .f = f_3d;

                    (
                        ghost_state(&j_state, end.gamma, normal, end.left),
                        Some(f_1d),
                    )
                }
                EndKind::Prescribed(bc) => {
                    let ghost = boundary_state(
                        &interior,
                        &end.reference,
                        bc,
                        end.area,
                        end.gamma,
                        end.left,
                    );

                    (ghost, Some(euler_flux(&ghost, end.gamma)))
                }
            };

            let pair = bc.get_mut(&end.pipe).expect("pipe has no boundary data");
            let side = match end.left {
                true => &mut pair.left,
                false => &mut pair.right,
            };
            side.ghost = Some(ghost);
            side.flux = flux;
        }
    }

    ///Advances every pipe and junction one full step, stage by stage.
    /// Leaves the decoded primitives stale, since q1 has moved.
    pub fn step(
        &mut self,
        dt: f64,
        pipes: &mut BTreeMap<usize, InteriorMethod>,
        junctions: &mut BTreeMap<usize, Junction>,
    ) {
        let stages = self.integrator.stages();

        //freeze q^n for everything before any stage runs
        for (id, pipe) in pipes.iter() {
            let regs = self.pipe_regs.get_mut(id).expect("pipe has no registers");
            regs.qn.copy_from(&pipe.solver().state().q1);
        }
        for (id, junction) in junctions.iter() {
            let regs = self
                .junction_regs
                .get_mut(id)
                .expect("junction has no registers");
            regs.qn = junction.q1;
        }

        for (k, &(a, b)) in stages.iter().enumerate() {
            let last = k + 1 == stages.len();

            self.refresh_boundary_data(k, junctions);
            self.stage_junctions(k, last, a, b, dt, junctions);
            self.stage_pipes(k, last, a, b, dt, pipes);
        }
    }

    ///Runs one stage for every junction.
    fn stage_junctions(
        &mut self,
        k: usize,
        last: bool,
        a: f64,
        b: f64,
        dt: f64,
        junctions: &mut BTreeMap<usize, Junction>,
    ) {
        for (id, junction) in junctions.iter_mut() {
            let regs = self
                .junction_regs
                .get_mut(id)
                .expect("junction has no registers");

            //evaluate the residual at the state this stage starts from
            let src = if k == 0 { regs.qn } else { regs.stages[k - 1] };
            junction.residual(&src);

            let next = blend(&regs.qn, &src, &junction.df, a, b, dt / junction.volume);
            if last {
                junction.q1 = next;
            } else {
                regs.stages[k] = next;
            }
        }
    }

    ///Runs one stage for every pipe.
    fn stage_pipes(
        &mut self,
        k: usize,
        last: bool,
        a: f64,
        b: f64,
        dt: f64,
        pipes: &mut BTreeMap<usize, InteriorMethod>,
    ) {
        for (id, pipe) in pipes.iter_mut() {
            //BoundaryPair is Copy, so take it before the registers are borrowed mutably
            let bc = *self.bc.get(id).expect("pipe has no boundary data");
            let regs = self.pipe_regs.get_mut(id).expect("pipe has no registers");

            //copy the scalars out first so the buffers below can be split-borrowed
            let s = pipe.solver().state();
            let (first, n_real, n_ghost) = (s.first, s.n_real, s.n_ghost);
            let c = dt / s.dx;

            //the ghosts belong to the state this stage reads, so they are filled from the
            //boundary data that was just resolved against it
            match k {
                0 => apply_bc(&mut regs.qn, n_ghost, &bc),
                _ => apply_bc(&mut regs.stages[k - 1], n_ghost, &bc),
            }

            //evaluate the residual at the state this stage starts from
            if k == 0 {
                pipe.solver_mut().residual(&regs.qn, &bc);
            } else {
                pipe.solver_mut().residual(&regs.stages[k - 1], &bc);
            }

            let Registers { qn, stages: bufs } = &mut *regs;

            if last {
                //last stage writes the answer straight into the pipe
                let PipeState { q1, df, .. } = pipe.solver_mut().state_mut();
                let src = if k == 0 { &*qn } else { &bufs[k - 1] };
                rk_stage(q1, qn, src, df, a, b, c, first, n_real);
                //so the stored state carries a boundary until the next step refreshes it
                apply_bc(q1, n_ghost, &bc);
            } else {
                let (done, rest) = bufs.split_at_mut(k);
                let dst = &mut rest[0];
                let src = if k == 0 { &*qn } else { &done[k - 1] };
                rk_stage(
                    dst,
                    qn,
                    src,
                    &pipe.solver().state().df,
                    a,
                    b,
                    c,
                    first,
                    n_real,
                );
            }
        }
    }
}

///Matches every pipe end marked as a junction to the inlet the junction holds for it.
fn resolve_ends(
    pipes: &BTreeMap<usize, InteriorMethod>,
    junctions: &BTreeMap<usize, Junction>,
) -> Vec<DrivenEnd> {
    let mut ends = Vec::new();

    for (&pipe, method) in pipes.iter() {
        let s = method.solver().state();

        for (bc, left) in [(s.left_bc, true), (s.right_bc, false)] {
            let kind = match bc {
                Some(BoundaryCondition::Junction(junction)) => {
                    let j = junctions
                        .get(&junction)
                        .expect("pipe names a junction that does not exist");
                    let index = j
                        .pipes
                        .iter()
                        .position(|inlet| inlet.pipe_id == pipe)
                        .expect("junction does not list the pipe that names it");
                    assert_eq!(
                        j.pipes[index].left, left,
                        "pipe {pipe} and junction {junction} disagree on which end they meet at"
                    );

                    EndKind::Junction {
                        junction,
                        inlet: index,
                        normal: j.pipes[index].normal,
                    }
                }
                Some(other) => EndKind::Prescribed(other),
                None => continue,
            };

            let cell = match left {
                true => s.first,
                false => s.first + s.n_real - 1,
            };

            ends.push(DrivenEnd {
                pipe,
                left,
                cell,
                reference: Vector3::new(s.q1[(0, cell)], s.q1[(1, cell)], s.q1[(2, cell)]),
                area: s.pipe_area(),
                gamma: s.gamma,
                kind,
            });
        }
    }

    ends
}
