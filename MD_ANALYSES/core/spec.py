#!/usr/bin/env python3
"""
spec.py -- what a system is, and which analyses it therefore deserves.

A spec is the one thing a person edits. It is proposed automatically from the
topology (core/detect.py) and confirmed by hand, in the same spirit as the
identifier and topic in the metadata pipeline: measured facts are filled in,
the judgements are yours.

    system:
      id: CDR/d2830a2a              # the record this belongs to
      label: 1F4                    # what to call it in a figure
      components: [glycoprotein, antibody, glycans, water, ions]
      bundle: antibody-complex      # which analyses follow from that
      selections:                   # VMD/CHARMM syntax; blank means the default
        receptor: segid A120 B120 C120 A41 B41 C41
        ligand:   segid IFH IFL IFM IFN IFQ IFS
        glycans:  segid GA1 to GA23 GAE GAH ...
        epitope:  segid A120 and resid 384 to 401
      mds:                          # one per replica, from the glossary
        - {id: 1-R1, topology: ..., trajectories: [...], framestep_ns: 0.2}

BUNDLES. A bundle is a list of stage names, nothing more, so supporting a new
kind of system is a table entry rather than new code. The components detected in
a topology select the bundle; the owner can override it.

STRIDES. Every stage declares the coarsest frame spacing at which its
observable still means something (`min_spacing_ns`). A stride that violates it
is refused with the spacing it would need -- at 0.2 ns/frame, stride 100 is one
frame per 20 ns, which is fine for occupancy and composition and useless for
exchange dynamics.
"""

import copy
import os

#: component -> the bundle it argues for, most specific last.
COMPONENT_BUNDLE = [
    ("protein", "apo-protein"),
    ("glycans", "glycoprotein"),
    ("nucleic", "nucleic-complex"),
    ("membrane", "membrane"),
    ("peptide", "peptide-complex"),
    ("small-molecule", "small-molecule"),
    ("antibody", "antibody-complex"),
]

#: bundle -> stages. `solvated` is added to any system that has water, and the
#: core set runs for everything.
BUNDLES = {
    "core": ["rmsd", "rmsf", "rg", "dssp", "pca"],
    "apo-protein": [],
    "glycoprotein": ["glycan_shielding", "glycan_hydration"],
    "antibody-complex": ["interface_contacts", "bsa", "key_residues",
                         "epitope_sasa", "binding_energy"],
    "small-molecule": ["ligand_rmsd", "binding_site", "contact_occupancy",
                       "binding_energy"],
    "peptide-complex": ["interface_contacts", "bsa", "key_residues",
                        "binding_energy"],
    "nucleic-complex": ["base_contacts", "groove_hydration",
                        "ion_condensation"],
    "membrane": ["leaflet", "insertion_depth"],
    "solvated": ["composition", "water_shell", "water_bridges",
                 "water_residence", "ions", "density_grids", "rdf",
                 "solvent_coupling"],
}

#: Which of those stage names actually exist as code today. The bundles above
#: are the PLAN; this is the part that runs. Keeping them apart matters: a
#: stage list that silently advertises `water_residence` when no such module
#: exists is a tool telling you it measured something it did not.
IMPLEMENTED = {"composition", "rmsd", "rmsf", "rg", "ions",
               "water_shell", "water_bridges", "water_residence",
               "density_grids", "rdf", "binding_energy",
               "interface_contacts"}


#: Default selections when a system declares none. Generic on purpose: they are
#: what "works on any system" means, and a system that needs more says so.
DEFAULT_SELECTIONS = {
    "solute": "not water and not ions",
    "protein": "protein",
    "glycans": "glycan",
    "water": "water",
    "ions": "ions",
    "heavy": "not name H*",
}

#: Stride defaults per stage family, in FRAMES.
DEFAULT_STRIDE = {"structure": 1, "surface": 5, "contacts": 5,
                  "solvent_static": 10, "solvent_dynamic": 1, "grids": 10}


class MD:
    """One replica: its topology, its segments in time order, its spacing."""

    __slots__ = ("id", "topology", "trajectories", "framestep_ns", "md_key",
                 "solvated_topology", "solvated_trajectories")

    def __init__(self, id, topology, trajectories, framestep_ns=None,
                 md_key=None, solvated_topology=None,
                 solvated_trajectories=None):
        self.id = id
        self.topology = topology
        self.trajectories = list(trajectories)
        self.framestep_ns = framestep_ns
        self.md_key = md_key or id
        # Structure stages read the stripped trajectory (55,818 atoms, 1.7 GB);
        # solvent stages must read the solvated one (767,697 atoms, 23 GB).
        # Both are kept, because no stage should have to guess which it needs.
        self.solvated_topology = solvated_topology
        self.solvated_trajectories = list(solvated_trajectories or [])

    def timeline_inputs(self):
        """The files whose HEADERS define this MD's time axis.

        Always the original segments, never a merged or stripped copy: cpptraj
        writes `merged.dcd` with istart=1, nsavc=1 and a placeholder delta, so
        a timeline built from it reports a 504 ns run as 0.003 ns and can only
        ever say "merged.dcd, frame 2016" -- which is true and useless. The
        merged file is a concatenation of the raw frames in order, so frame i
        of it IS raw frame i of these segments, and this mapping is what lets a
        click on a plot name `pr.part6.dcd, local frame 251`.
        """
        return self.solvated_trajectories or self.trajectories

    def inputs(self, solvated=False):
        if solvated:
            if not self.solvated_trajectories:
                raise ValueError(
                    f"{self.id}: this analysis needs the SOLVATED trajectory, "
                    f"and none is declared. Water and ions were stripped from "
                    f"{os.path.basename(self.trajectories[0])}.")
            return self.solvated_topology or self.topology, \
                self.solvated_trajectories
        return self.topology, self.trajectories

    def as_dict(self):
        return {"id": self.id, "md_key": self.md_key,
                "topology": self.topology, "trajectories": self.trajectories,
                "framestep_ns": self.framestep_ns,
                "solvated_topology": self.solvated_topology,
                "solvated_trajectories": self.solvated_trajectories}


class Spec:
    """One system: what it is made of, how to select its parts, and its MDs."""

    def __init__(self, id, label=None, components=(), bundle=None,
                 selections=None, mds=(), stride=None, options=None):
        self.id = id
        self.label = label or id.split("/")[-1]
        self.components = list(components)
        self.bundle = bundle or bundle_for(self.components)
        self.selections = dict(DEFAULT_SELECTIONS)
        supplied = {k: v for k, v in (selections or {}).items()
                    if str(v or "").strip()}
        self.selections.update(supplied)
        #: The names a PERSON wrote, as opposed to the generic defaults above.
        #: A stage needs the difference: `glycans: glycan` matching nothing
        #: means this system has no glycans, which is ordinary; the same
        #: selection matching nothing after someone typed `segid GA1 to GA23`
        #: means they typed it wrong, which is not.
        self.declared = set(supplied)
        self.mds = list(mds)
        self.stride = dict(DEFAULT_STRIDE)
        self.stride.update(stride or {})
        self.options = dict(options or {})

    # ------------------------------------------------------------------ stages
    def stages(self, implemented_only=False):
        """Every stage this system's components ask for, in a stable order.

        `implemented_only` gives the subset that exists as code; the rest is
        the roadmap, and `planned()` names it so the difference is visible
        rather than implied.
        """
        names = list(BUNDLES["core"])
        for b in (self.bundle, "solvated" if self.has_solvent else None):
            for s in BUNDLES.get(b or "", []):
                if s not in names:
                    names.append(s)
        if implemented_only:
            return [n for n in names if n in IMPLEMENTED]
        return names

    def planned(self):
        """Stages this system would get, once they are written."""
        return [n for n in self.stages() if n not in IMPLEMENTED]

    @property
    def has_solvent(self):
        return "water" in self.components or "ions" in self.components

    def selection(self, name):
        return self.selections.get(name, DEFAULT_SELECTIONS.get(name, ""))

    def stride_for(self, family):
        return int(self.stride.get(family, 1))

    # -------------------------------------------------------------------- i/o
    def as_dict(self):
        return {"id": self.id, "label": self.label,
                "components": self.components, "bundle": self.bundle,
                "selections": {k: v for k, v in self.selections.items()},
                "declared": sorted(self.declared),
                "stride": self.stride, "options": self.options,
                "mds": [m.as_dict() for m in self.mds]}

    @classmethod
    def from_dict(cls, d):
        mds = [MD(**{k: v for k, v in m.items() if k in MD.__slots__})
               for m in d.get("mds", [])]
        spec = cls(id=d["id"], label=d.get("label"),
                   components=d.get("components", []),
                   bundle=d.get("bundle"), selections=d.get("selections"),
                   mds=mds, stride=d.get("stride"), options=d.get("options"))
        # A round trip hands every selection back as if it had been supplied,
        # so the declared set is restored from the record rather than inferred.
        if "declared" in d:
            spec.declared = set(d["declared"] or ())
        return spec

    def copy(self):
        return Spec.from_dict(copy.deepcopy(self.as_dict()))

    def __repr__(self):                                    # pragma: no cover
        return (f"<Spec {self.id} {self.bundle} "
                f"{len(self.mds)} MD(s) {self.components}>")


def bundle_for(components):
    """The most specific bundle the components argue for."""
    chosen = "apo-protein"
    for comp, bundle in COMPONENT_BUNDLE:
        if comp in components:
            chosen = bundle
    return chosen


class SpacingTooCoarse(ValueError):
    """An analysis was asked for at a spacing where it means nothing."""


def check_spacing(stage, framestep_ns, stride, min_spacing_ns):
    """Refuse an observable the sampling cannot support, and say what it needs.

    Residence times, diffusion and hydrogen-bond lifetimes at one frame per
    20 ns are not approximate -- they are meaningless, and they look perfectly
    plausible in a table. The number this returns is the effective spacing.
    """
    if not framestep_ns:
        return None
    effective = float(framestep_ns) * max(1, int(stride))
    if min_spacing_ns and effective > float(min_spacing_ns):
        raise SpacingTooCoarse(
            f"{stage}: frames are {effective:g} ns apart "
            f"({framestep_ns:g} ns x stride {stride}), but this observable "
            f"needs {min_spacing_ns:g} ns or finer. Lower the stride, or run "
            f"a window with denser output.")
    return effective
