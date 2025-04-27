# Dune_NonSaturated
Finite element model to simulate problems involving the coupling between a fluid and a non-cohesive layer: application to the migration of a barchan dune.

Inputs:

  - MESH.unv: Finite element mesh with boundary conditions, in .unv format.

  - alpha_t0.dat: File containing parameters related to the fluid, the sediment, and the simulation settings.

Outputs:

  - nsvtk2D_*.vtk**: Fluid solutions in ParaView format, including the pressure and velocity fields, among others.

  - nsvtk2D_Bed_*.vtk**: Fluid-sediment interface location.

  - nsvtk2D_BedUndef_*.vtk**: Contains data of the saltation layer, such as the sediment flux, grain velocity, and saltation layer density. The fluid-sediment interface is represented as a flat line; these files are suitable for plotting saltation layer data.
