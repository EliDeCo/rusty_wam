use crate::boundaries::{BoundaryCondition, boundary_state, euler_flux};
use crate::junctions::{Junction, ghost_state, interface_flux};
use crate::pipes::{BoundaryPair, InteriorMethod, PipeState, TimeIntegrator, apply_bc, rk_stage};
use nalgebra::{Matrix3xX, Vector3, Vector5};
use std::collections::BTreeMap;

///Scratch buffers one pipe needs to be walked through the integrator's stages.
struct Registers {
    ///id of the pipe these belong to, which only a debug assertion reads
    id: usize,
    qn: Matrix3xX<f64>,
    stages: Vec<Matrix3xX<f64>>,
}

///The same buffers for a junction, which is a single cell of five values.
struct JunctionRegisters {
    ///id of the junction these belong to, which only a debug assertion reads
    id: usize,
    qn: Vector5<f64>,
    stages: Vec<Vector5<f64>>,
}

///One integrator stage for a junction cell: dst = a*q_n + b*src - c*dt/vol*df.
/// The junction is a control volume, so its scaling is dt/volume rather than dt/dx.
#[allow(clippy::too_many_arguments)]
fn blend(
    q_n: &Vector5<f64>,
    src: &Vector5<f64>,
    df: &Vector5<f64>,
    a: f64,
    b: f64,
    c: f64,
    dt_over_vol: f64,
) -> Vector5<f64> {
    debug_assert!((a + b - 1.0).abs() < 1e-15, "stage weights must sum to 1");

    a * q_n + b * src - (c * dt_over_vol) * df
}

///What drives one pipe end, resolved once so no stage has to search for it.
enum EndKind {
    Junction {
        ///id of the junction, which is how the inlet this flux belongs to is reached
        junction: usize,
        ///slot the same junction holds in the driver's register list
        slot: usize,
        ///index of this pipe's entry in the junction's own inlet list
        inlet: usize,
        normal: Vector3<f64>,
    },
    Prescribed(BoundaryCondition),
}

///One pipe end whose ghosts the driver fills every stage.
struct DrivenEnd {
    ///slot the pipe holds in the register and boundary lists
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
    //held in the order the pipes and junctions iterate, which is fixed once Driver::new
    //has run, so a stage indexes straight in rather than searching by id
    pipe_regs: Vec<Registers>,
    junction_regs: Vec<JunctionRegisters>,
    ///what each junction hands its branches, rebuilt every stage
    bc: Vec<BoundaryPair>,
    ends: Vec<DrivenEnd>,
}

impl Driver {
    ///Allocates every stage register once, pairing each pipe's method with the integrator
    /// it is cheapest under. One stage loop drives every pipe and junction, so the network
    /// takes a single integrator: the most capable if the methods ever disagree, since a
    /// method that wants Euler is stable under Ssp3 and never the reverse.
    pub fn new(
        pipes: &BTreeMap<usize, InteriorMethod>,
        junctions: &BTreeMap<usize, Junction>,
    ) -> Self {
        let integrator = pipes
            .values()
            .map(|p| p.kind().integrator())
            .max_by_key(|i| i.stages().len())
            .unwrap_or(TimeIntegrator::Euler);

        Self::with_integrator(integrator, pipes, junctions)
    }

    ///The same, with the integrator named rather than derived. Only the scratch harnesses
    /// use this, to sweep pairings the shipped solver would never select.
    pub fn with_integrator(
        integrator: TimeIntegrator,
        pipes: &BTreeMap<usize, InteriorMethod>,
        junctions: &BTreeMap<usize, Junction>,
    ) -> Self {
        //forward Euler is unstable with any reconstruction at every Courant number, at
        //kappa = 0 as well as 1/3, and fails thousands of steps in rather than at once
        let reconstructs = pipes.values().any(|p| {
            matches!(
                p,
                InteriorMethod::Muscl2RoeM1D(_) | InteriorMethod::Muscl3RoeM1D(_)
            )
        });
        assert!(
            !(matches!(integrator, TimeIntegrator::Euler) && reconstructs),
            "a MUSCL method is unstable with Euler at every Courant number; use Ssp2 or higher"
        );

        let n_regs = integrator.n_registers();

        let pipe_regs = pipes
            .iter()
            .map(|(&id, pipe)| {
                let n_total = pipe.solver().state().n_total;

                Registers {
                    id,
                    qn: Matrix3xX::zeros(n_total),
                    stages: (0..n_regs).map(|_| Matrix3xX::zeros(n_total)).collect(),
                }
            })
            .collect();

        let junction_regs = junctions
            .keys()
            .map(|&id| JunctionRegisters {
                id,
                qn: Vector5::zeros(),
                stages: vec![Vector5::zeros(); n_regs],
            })
            .collect();

        let bc = vec![BoundaryPair::default(); pipes.len()];
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
        for pair in self.bc.iter_mut() {
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
            let p_regs = &pipe_regs[end.pipe];
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
                    slot,
                    inlet,
                    normal,
                } => {
                    let j_regs = &junction_regs[*slot];
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

            let pair = &mut bc[end.pipe];
            let side = match end.left {
                true => &mut pair.left,
                false => &mut pair.right,
            };
            side.ghost = Some(ghost);
            side.flux = flux;
        }
    }

    ///Freezes the state the coming step starts from, resolves every boundary against it
    /// and fills the ghosts, so `get_timestep` decodes what the first stage will read.
    pub fn prepare(
        &mut self,
        pipes: &mut BTreeMap<usize, InteriorMethod>,
        junctions: &mut BTreeMap<usize, Junction>,
    ) {
        //freeze q^n for everything before any stage runs
        for (regs, (&id, pipe)) in self.pipe_regs.iter_mut().zip(pipes.iter()) {
            debug_assert_eq!(regs.id, id, "pipe registers are out of order");
            regs.qn.copy_from(&pipe.solver().state().q1);
        }
        for (regs, junction) in self.junction_regs.iter_mut().zip(junctions.values()) {
            debug_assert_eq!(regs.id, junction.id, "junction registers are out of order");
            regs.qn = junction.q1;
        }

        self.refresh_boundary_data(0, junctions);

        //q1 gets the same ghosts as the register it was frozen into, which is what makes
        //the decode behind get_timestep serve the first stage as well
        let Self { pipe_regs, bc, .. } = self;
        for ((regs, bc), pipe) in pipe_regs.iter_mut().zip(bc.iter()).zip(pipes.values_mut()) {
            let n_ghost = pipe.solver().state().n_ghost;
            apply_bc(&mut pipe.solver_mut().state_mut().q1, n_ghost, bc);
            apply_bc(&mut regs.qn, n_ghost, bc);
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

        for (k, &(a, b, c)) in stages.iter().enumerate() {
            let last = k + 1 == stages.len();

            //stage 0 reads what prepare resolved, so only later stages need refreshing
            if k > 0 {
                self.refresh_boundary_data(k, junctions);
            }
            self.stage_junctions(k, last, a, b, c, dt, junctions);
            self.stage_pipes(k, last, a, b, c, dt, pipes);
        }
    }

    ///Runs one stage for every junction.
    #[allow(clippy::too_many_arguments)]
    fn stage_junctions(
        &mut self,
        k: usize,
        last: bool,
        a: f64,
        b: f64,
        c: f64,
        dt: f64,
        junctions: &mut BTreeMap<usize, Junction>,
    ) {
        for (regs, junction) in self.junction_regs.iter_mut().zip(junctions.values_mut()) {
            debug_assert_eq!(regs.id, junction.id, "junction registers are out of order");

            //evaluate the residual at the state this stage starts from, which at stage 0
            //is the state get_timestep already decoded
            let src = if k == 0 { regs.qn } else { regs.stages[k - 1] };
            if k > 0 {
                junction.decode_from(&src);
            }
            junction.residual(&src);

            let next = blend(&regs.qn, &src, &junction.df, a, b, c, dt / junction.volume);
            if last {
                junction.q1 = next;
            } else {
                regs.stages[k] = next;
            }
        }
    }

    ///Runs one stage for every pipe.
    #[allow(clippy::too_many_arguments)]
    fn stage_pipes(
        &mut self,
        k: usize,
        last: bool,
        a: f64,
        b: f64,
        c: f64,
        dt: f64,
        pipes: &mut BTreeMap<usize, InteriorMethod>,
    ) {
        let Self { pipe_regs, bc, .. } = self;
        for ((regs, bc), pipe) in pipe_regs.iter_mut().zip(bc.iter()).zip(pipes.values_mut()) {
            //copy the scalars out first so the buffers below can be split-borrowed
            let s = pipe.solver().state();
            let (first, n_real, n_ghost) = (s.first, s.n_real, s.n_ghost);
            let dt_over_dx = dt / s.dx;

            //the ghosts belong to the state this stage reads, so they are filled from the
            //boundary data that was just resolved against it
            if k > 0 {
                apply_bc(&mut regs.stages[k - 1], n_ghost, bc);
            }

            //evaluate the residual at the state this stage starts from, whose primitives
            //stage 0 already has from get_timestep
            if k == 0 {
                pipe.solver_mut().residual(&regs.qn, bc, true);
            } else {
                pipe.solver_mut().residual(&regs.stages[k - 1], bc, false);
            }

            let Registers {
                qn, stages: bufs, ..
            } = &mut *regs;

            if last {
                //last stage writes the answer straight into the pipe
                let PipeState { q1, df, .. } = pipe.solver_mut().state_mut();
                let src = if k == 0 { &*qn } else { &bufs[k - 1] };
                rk_stage(q1, qn, src, df, a, b, c, dt_over_dx, first, n_real);
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
                    dt_over_dx,
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
    //the slot each junction holds in the driver's register list, which is this map's order
    let slots: BTreeMap<usize, usize> = junctions
        .keys()
        .enumerate()
        .map(|(slot, &id)| (id, slot))
        .collect();

    let mut ends = Vec::new();

    for (pipe_slot, (&pipe, method)) in pipes.iter().enumerate() {
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
                        slot: slots[&junction],
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
                pipe: pipe_slot,
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
