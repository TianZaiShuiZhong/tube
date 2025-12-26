from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Node:
    id: int
    x: float
    y: float
    z: float


@dataclass(frozen=True)
class Element:
    id: int
    etype: int  # 288/290 etc
    n1: int
    n2: int
    n3: int | None = None
    n4: int | None = None
    mat: int | None = None
    sec: int | None = None

    def end_nodes(self) -> tuple[int, int]:
        return (self.n1, self.n2)

    def orientation_node(self) -> int | None:
        return self.n3


@dataclass(frozen=True)
class SectionPipe:
    id: int
    outer_diameter: float
    thickness: float

    @property
    def ro(self) -> float:
        return 0.5 * self.outer_diameter

    @property
    def ri(self) -> float:
        return self.ro - self.thickness


@dataclass(frozen=True)
class PipeModalConfig:
    # Low-order circumferential Fourier modes: m=0 (breathing/axisym), m=1 (bending), m=2 (oval)
    modes: tuple[int, ...] = (0,)  # typically (0,1,2) for full capability
    # number of circumferential divisions for visualization/expansion (default 24 or 36)
    n_circ_vis: int = 24
    # amplitude scaling for modal visualization (default 1.0)
    modal_scale: float = 1.0


@dataclass(frozen=True)
class ModalDOF:
    # A single modal DOF entry for a node and circumferential mode
    node: int
    mode: int  # Fourier index (0=axisymmetric, 1=bending, 2=oval, etc.)
    # DOF type identifiers (examples: 'ur','uz','ut','warp','ellip') to be refined
    dof: str


@dataclass(frozen=True)
class FourierModalState:
    # Per-node modal amplitude state for visualization/output
    # Format: {node_id: {mode: amplitude}, ...}
    # Example: {1: {0: 0.5, 1: 0.1, 2: 0.0}, 2: {0: 0.5, 1: -0.1, 2: 0.05}, ...}
    modal_amplitudes: dict[int, dict[int, float]] = field(default_factory=dict)


@dataclass(frozen=True)
class Material:
    id: int
    ex: float
    nuxy: float
    dens: float | None = None
    alp: float | None = None
    reft: float | None = None
    creep: tuple[float, float, float, float] | None = None
    plasticity: tuple[float, float] | None = None  # (yield_strength, tangent_modulus)


@dataclass(frozen=True)
class Constraint:
    node: int
    dof: str
    value: float


@dataclass(frozen=True)
class NodalLoad:
    node: int
    dof: str
    value: float


@dataclass(frozen=True)
class ElementPressure:
    elem: int
    value: float


@dataclass
class Model:
    nodes: dict[int, Node] = field(default_factory=dict)
    elements: dict[int, Element] = field(default_factory=dict)
    materials: dict[int, Material] = field(default_factory=dict)
    sections: dict[int, SectionPipe] = field(default_factory=dict)

    constraints: list[Constraint] = field(default_factory=list)
    nodal_loads: list[NodalLoad] = field(default_factory=list)
    elem_pressures: list[ElementPressure] = field(default_factory=list)

    uniform_temperature: float | None = None
    accel: tuple[float, float, float] | None = None

    # analysis controls (from CDB)
    time_end: float | None = None
    nsubsteps: int | None = None

    # modal configuration (scaffolding)
    modal_config: PipeModalConfig = field(default_factory=PipeModalConfig)
    modal_dofs: list[ModalDOF] = field(default_factory=list)
    # modal state for visualization (populated during solve or postproc)
    modal_state: FourierModalState = field(default_factory=FourierModalState)
