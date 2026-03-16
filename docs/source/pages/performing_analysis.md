# Loads, load cases, and analysis

*ospgrillage* contains a load module which wraps `OpenSeesPy` commands to define loads, assemble load cases, and perform analysis.

For all example code in this page, *ospgrillage* is imported as `og`

```python
import ospgrillage as og
```

## Load analysis workflow

Figure 1 shows the flowchart for the load module of *ospgrillage*.

![Figure 1: Load analysis utility flow chart](../images/analysis_workflow.png)

## Defining loads

Loads are created with the interface function {func}`~ospgrillage.load.create_load`. Users pass the `loadtype` keyword argument to specify the load type. Available load types include {ref}`point`, {ref}`line`, and {ref}`patch` loads.

Each load type requires the user to specify its load point(s). This is achieved by the {func}`~ospgrillage.load.create_load_vertex` function. This function creates a {class}`~ospgrillage.load.LoadPoint` `namedtuple` `(x, y, z, p)` where `x`, `y`, `z` are the coordinates of the load point and `p` is the magnitude of the vertical loading. The interpretation of `p` depends on the load type (explained below). By default, `y` is `0` — i.e. the grillage model plane.

```python
point_load_location = og.create_load_vertex(x=5, z=2, p=20)  # create load point
```

Depending on the load type, a minimum number of {class}`~ospgrillage.load.LoadPoint` tuples are required. These are passed to the `point1`, `point2`, … keyword arguments of {func}`~ospgrillage.load.create_load`. The following sections explain the required vertices for each load type.

Loads are generally defined in the global coordinate system with respect to the created grillage model. However, a user-defined local coordinate system is required when defining {ref}`compound-load` later on.

### Nodal loads

Nodal loads are loads applied directly onto nodes of the grillage model. Nodal loads are defined using {func}`~ospgrillage.load.create_load`, specifying `loadtype="nodal"`. There are six degrees-of-freedom (DOFs) for acting loads at each node.

Nodal loads do not require a load vertex. Instead, they require a `NodeForces(Fx, Fy, Fz, Mx, My, Mz)` `namedtuple`. The following example creates a `NodeForces` tuple and a nodal load on Node 13 of a model, with 10 unit force in both the X and Y directions.

```python
nodalforce = og.NodeForces(Fx=10, Fy=10)
node13force = og.create_load(loadtype="nodal", name="nodal 13", node_force=nodalforce)  # other DOFs default to 0
```

```{note}
You only need to specify non-zero values for the desired DOFs. Any unspecified component defaults to zero.
```

(point)=
### Point Loads

A point load is a force applied at a single infinitesimal point of the grillage model. Point loads can represent a wide range of loads such as truck axles or superimposed dead loads on a deck.

Point loads are created using {func}`~ospgrillage.load.create_load`, passing `loadtype="point"`. A point load takes only a single {class}`~ospgrillage.load.LoadPoint` tuple. `p` should have units of force (e.g. N, kN, kips) — see Figure 2.

![Figure 2: Point load](../images/point.png)

The following example code creates a 20 force unit point load located at (5,0,2) in the global coordinate system.

```python
point_load_location = og.create_load_vertex(x=5, z=2, p=20)  # create load point
point_load = og.create_load(loadtype="point",name="single point", point1=point_load_location)
```

(line)=
### Line Loads

Line loads are loads exerted along a line. They are useful to represent loads such as self-weight of longitudinal beams or distributed loads along beam elements.

Line loads are created with {func}`~ospgrillage.load.create_load` passing `loadtype="line"` and require at least two {class}`~ospgrillage.load.LoadPoint` values (corresponding to the start and end of the line) — see Figure 3. Using more than two points allows a varying line-load profile. `p` should have units of force per distance (e.g. kN/m, kips/ft).

![Figure 3: Line load](../images/line.png)

The following example code is a constant Two force per distance unit line load (UDL) in the global coordinate system from -1 to 11 distance units in the `x`-axis and along the position in the `z`-axis at 3 distance units.

```python
barrier_point_1 = og.create_load_vertex(x=-1, z=3, p=2)
barrier_point_2 = og.create_load_vertex(x=11, z=3, p=2)
Barrier = og.create_load(loadtype="line", name="Barrier curb", point1=barrier_point_1, point2=barrier_point_2)
```

```{note}
As of release 0.1.0, curved line loads are not available.
```

(patch)=
### Patch loads

Patch loads represent loads distributed uniformly over an area, such as traffic lanes.

Patch loads are created with {func}`~ospgrillage.load.create_load`, specifying `loadtype="patch"`. A patch load requires at least four {class}`~ospgrillage.load.LoadPoint` tuples (corresponding to the vertices of the loaded area) — see Figure 4. Using eight tuples allows a curved surface loading profile. `p` should have units of force per area (e.g. kN/m², kips/ft²).

![Figure 4: Patch load](../images/patch.png)

The following example code creates a constant 5 force per area unit patch load in the global coordinate system.

```python
lane_point_1 = og.create_load_vertex(x=0, z=3, p=5)
lane_point_2 = og.create_load_vertex(x=8, z=3, p=5)
lane_point_3 = og.create_load_vertex(x=8, z=5, p=5)
lane_point_4 = og.create_load_vertex(x=0, z=5, p=5)
Lane = og.create_load(loadtype="patch",name="Lane 1", point1=lane_point_1, point2=lane_point_2, point3=lane_point_3, point4=lane_point_4)
```

```{note}
As of release 0.1.0, curved patch loads are not available.
```

(compound-load)=
## Compound loads

Two or more of the basic load types can be combined to form a Compound load. All load types are applied in the direction of the global $y$-axis. Loads in other directions and applied moments are currently not supported.

To create a compound load, use the {func}`~ospgrillage.load.create_compound_load` function. This function creates a {class}`~ospgrillage.load.CompoundLoad` object.

Compound load are typically defined in a **local coordinate system** and then set to global coordinate system of the grillage. Figure 5 shows the relationship and process of mapping local to global system of a compound load.

![Figure 5: Compound load](../images/compoundload.png)

The following code creates a point and line load which is to be assigned as a Compound load.

```python
# components in a compound load
wheel_1 = og.create_load(loadtype="point", point1= og.create_load_vertex(x=0, z=3, p=5))  # point load 1
wheel_2 = og.create_load(loadtype="point", point1= og.create_load_vertex(x=0, z=3, p=5))  # point load 2
```

The following code creates a Compound load and adds the created load objects to it:

```python
C_Load = og.create_compound_load(name = "Axle tandem")  # constructor of compound load
C_Load.add_load(load_obj=wheel_1) # add wheel_1
C_Load.add_load(load_obj=wheel_2) # add wheel_2
```

After defining all required load objects, {class}`~ospgrillage.load.CompoundLoad` requires users to define the global coordinate to map the origin of the local coordinate system to the global coordinate space. This is done using {meth}`~ospgrillage.load.CompoundLoad.set_global_coord`, passing a {class}`~ospgrillage.mesh.Point` `namedtuple` `(x, y, z)` as seen in Figure 5. If not specified, the mapping's reference point defaults to the **Origin** (0, 0, 0).

The following example sets the local **Origin** of the compound load, including all load points for all load objects of **C_load** by x + 4, y + 0, and z + 3.

```python
C_Load.set_global_coord(Point(4,0,3))
```

**Coordinate System**

When adding each load object, the {class}`~ospgrillage.load.CompoundLoad` class allow users to input a `load_coord=` keyword argument. This relates to the load object - whether it was previously defined in the user-defined *local* or in the *global* coordinate system. The following explains the various input conditions

```{note}
Compound loads require users to pay attention to the difference between local and global coordinate systems (see {doc}`package_design` for more information).

The {class}`~ospgrillage.load.CompoundLoad` stores load objects in the **local coordinate system**. When passed to a {class}`~ospgrillage.load.LoadCase`, the compound load's vertices are automatically converted to **global coordinates** based on the {meth}`~ospgrillage.load.CompoundLoad.set_global_coord` mapping.
```

(load-cases)=
## Load cases

Load cases are a set of load types ({ref}`point`, {ref}`line`, {ref}`patch`, {ref}`compound-load`) used to define a particular loading condition. Compound loads are treated as a single load group within a load case having same reference points (e.g. tandem axle) and properties (e.g. load factor)

After load objects are created, users add them to {class}`~ospgrillage.load.LoadCase` objects. First, instantiate a {class}`~ospgrillage.load.LoadCase` using {func}`~ospgrillage.load.create_load_case`:

```python
DL = og.create_load_case(name="Dead Load")
```

Then add load objects using {meth}`~ospgrillage.load.LoadCase.add_load`. The following code shows how the above load types are added to the *DL* load case.

```python
DL.add_load(point_load)  # each line adds individual load types to the load case
DL.add_load(Barrier)
DL.add_load(Lane)
```

After adding loads, the {class}`~ospgrillage.load.LoadCase` object is added to the grillage model for analysis using {meth}`~ospgrillage.osp_grillage.OspGrillage.add_load_case`. Repeat this step for each defined load case.

```python
example_bridge.add_load_case(DL)  # adding this load case to grillage model
```

(moving-load)=
## Moving load

For moving load analysis, users create moving load objects using {class}`~ospgrillage.load.MovingLoad` class. The moving load class takes a load type object ({ref}`point`, {ref}`line`, {ref}`patch`, {ref}`compound-load`) and moves the load through a path points described by a {class}`~ospgrillage.load.Path` object.

Figure 6 summarizes the relationship between moving loads, paths and the position of the loads on the grillage model.

![Figure 6: Moving load](../images/movingload.png)

### Moving path

{class}`~ospgrillage.load.Path` object is created using {func}`~ospgrillage.load.create_moving_path`.

{class}`~ospgrillage.load.Path` requires two {class}`~ospgrillage.mesh.Point` `namedtuples` `(x, y, z)` to describe its start and end position. The following example creates a path from 2 to 4 distance units in the global coordinate system.

```python
single_path = og.create_moving_path(start_point=og.Point(2,0,2), end_point= og.Point(4,0,2))
```

### Creating moving load

The following example code creates a compound load consisting of two point loads moving along the defined **single_path**

```python
# create components of compound load
front_wheel = og.create_load_vertex(x=0, z=0, p=6)
back_wheel = og.create_load_vertex(x=-1, z=0, p=6)
Line = og.create_load(loadtype="line",point1=front_wheel,point2=back_wheel)
tandem = og.create_compound_load("Two wheel vehicle")

move_line = og.create_moving_load(name="Line Load moving") # moving load obj
move_line.set_path(single_path)   # set path
move_line.add_load(load_obj=Line)  # add compound load to moving load
```

From here, use {meth}`~ospgrillage.osp_grillage.OspGrillage.add_load_case` to add the moving load to the grillage model. The function automatically creates multiple incremental {ref}`load-cases`, each corresponding to a load position along the moving path.

```python
example_bridge.add_load_case(move_point)
```

### Advanced usage

All basic loads added to a {class}`~ospgrillage.load.MovingLoad` via {meth}`~ospgrillage.load.MovingLoad.add_load` are assigned a single common {class}`~ospgrillage.load.Path`.

{class}`~ospgrillage.load.MovingLoad` also supports individual paths for each load. For this, pass a `global_increment` parameter when creating the moving load to ensure all paths use the same increment. Then, pass the `path_obj` keyword argument when adding each load.

The following example shows this procedure:

```python
# create moving load with global increment of 20 for all unique moving path
moving_load_group = og.create_moving_load(name="Line Load moving",global_increment=20)

# add load + their respective path
move_load_group.add_load(load_obj=truck_a,path_obj=path_a)
move_load_group.add_load(load_obj=truck_b,path_obj=path_b)
```

## Running analysis

Once all defined load cases (static and moving) have been added to the grillage object, analysis can be conducted.

To analyse load case(s), call {meth}`~ospgrillage.osp_grillage.OspGrillage.analyze`. By default this runs all defined load cases. To run only specific load cases, pass a load case name `str` or a `list` of names to the `load_case` keyword argument. The following example shows the various options:

```python
# analyze all
example_bridge.analyze()
# or a single str
example_bridge.analyze(load_case="DL")
# or a single element list
example_bridge.analyze(load_case=["DL"])
# or a list of multiple load cases
example_bridge.analyze(load_case=["DL","SDL"])
```
